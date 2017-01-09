#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>

#include <fftw3.h>
#include <x86intrin.h>

static int packet_fd;
static int read_fd;

static bool enable_frequency_tracking = true;
static bool capture_qpsk = false;
static bool capture_fft = false;

static unsigned int interleave_table[1536];
static uint32_t depunctured_data[384*4];

static unsigned int depuncture(unsigned int start, unsigned int stop,
	uint32_t depuncture_word, uint32_t *data_in, uint32_t *data_out)
{
	uint64_t word_in;

#if 1
	unsigned int len;
	unsigned int level;
	unsigned int num_words;

	data_in = data_in + start / 32;
	word_in = *data_in++;
	word_in |= ((uint64_t)(*data_in++)) << 32;
	word_in >>= (start % 32);
	level = 64 - (start % 32);

	len = __builtin_popcount(depuncture_word);

	num_words = (stop - start) / len;

	while (start < stop) {
		*data_out++ = __builtin_ia32_pdep_si((uint32_t)word_in, depuncture_word);
		word_in >>= len;
		start += len;
		level -= len;
		if (level < 32) {
			word_in |= ((uint64_t)(*data_in++)) << level;
			level += 32;
		}
	}

	return num_words;
#else
	uint32_t word_out = 0;
	unsigned int n, o, m;
	uint32_t x = 0x1;
	uint32_t bit;

	n = 0;
	m = 0;

	while (start < stop || n != 0) {
		if (depuncture_word & x) {
			o = start & 0x1f;
			bit = (word_in >> o) & 1;
			word_out |= bit << n;
			start++;
			if (o == 31) {
				data_in++;
				word_in = *data_in;
			}
		}

		n++;
		if (n == 32) {
			x = 0x1;
			*data_out++ = word_out;
			word_out = 0;
			n = 0;
			m++;
		} else {
			x <<= 1;
		}
	}
//	printf("%d\n", m);
#endif
}

struct viterbi_state {
	unsigned int num_paths;
	unsigned int path_metric_index;
	uint16_t path_metric[64];
	uint64_t path_choices[51200];
};

struct viterbi_state vb;

static void viterbi_reset()
{
	vb.num_paths = 0;
	vb.path_metric_index = 0;
	memset(vb.path_metric, 0x0, sizeof(vb.path_metric));
}

static __attribute__ ((noinline))  void viterbi_calc_branch_metric(uint32_t data, uint32_t puncture_word,
	uint16_t (*branch_metric)[8])
{
	unsigned int i, j;
	uint32_t x, y;

	for (i = 0; i < 8; i++) {
		x = (i << 1) | (i >> 2);
		x |= x << 4;
		x |= x << 8;
		x |= x << 16;
		y = (x ^ data) & puncture_word;
		y = (y & 0x55555555) + ((y & 0xaaaaaaaa) >> 1);
		y = (y & 0x33333333) + ((y & 0xcccccccc) >> 2);
		for (j = 0; j < 8; j++)
			branch_metric[j][i] = (y >> (4*j)) & 0x3;
	}
}

typedef uint16_t v8hu __attribute__((vector_size(16)));
typedef int16_t v8hi __attribute__((vector_size(16)));

static __attribute__ ((noinline)) v8hu vpmin(v8hu va, v8hu vb)
{
#if 0
	return (v8hu)__builtin_ia32_pminuw128((v8hi)va, (v8hi)vb);
#else
	uint16_t *a, *b;
	uint16_t c[8];
	unsigned int i;

	a = (uint16_t *)&va;
	b = (uint16_t *)&vb;

	for (i = 0; i < 8; i++)
		c[i] = (a[i] < b[i]) ? a[i] : b[i];
	return *(v8hu*)&c;
#endif
}

static void update_path_metric(v8hu *bm,
	v8hu *path_metric, v8hu *path_metric_next,
	unsigned int k, uint64_t *pc)
{
	v8hu vbm_shuffled;
	v8hu vpm_a, vpm_b;
	v8hu compare0;
	v8hu compare1;
	v8hu vbm;
	uint16_t path_choice = 0;
	unsigned int i;

	static const v8hu shuffle_mask0 = {0, 0, 1, 1, 2, 2, 3, 3};
	static const v8hu shuffle_mask1 = {4, 4, 5, 5, 6, 6, 7, 7};
	static const v8hu shuffle_mask_swap = {1, 0, 3, 2, 5, 4, 7, 6};

	static const v8hu shuffle_mask_bm[8] = {
		{0, 7, 3, 4, 5, 2, 6, 1},
		{5, 2, 6, 1, 0, 7, 3, 4},
		{2, 5, 1, 6, 7, 0, 4, 3},
		{7, 0, 4, 3, 2, 5, 1, 6},
		{4, 3, 7, 0, 1, 6, 2, 5},
		{1, 6, 2, 5, 4, 3, 7, 0},
		{6, 1, 5, 2, 3, 4, 0, 7},
		{3, 4, 0, 7, 6, 1, 5, 2}
	};

	vbm = *bm;

	vbm_shuffled = __builtin_shuffle(vbm, shuffle_mask_bm[2*k]);
	vpm_a = __builtin_shuffle(path_metric[k], shuffle_mask0);
	vpm_a = vpm_a + vbm_shuffled;

	vbm_shuffled = __builtin_shuffle(vbm_shuffled, shuffle_mask_swap);
	vpm_b = __builtin_shuffle(path_metric[4+k], shuffle_mask0);
	vpm_b = vpm_b + vbm_shuffled;

	path_metric_next[2*k] = vpmin(vpm_a, vpm_b);

	compare0 = vpm_b < vpm_a;
	compare0 &= (v8hu){0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};

	vbm_shuffled = __builtin_shuffle(vbm, shuffle_mask_bm[2*k+1]);
	vpm_a = __builtin_shuffle(path_metric[k], shuffle_mask1);
	vpm_a = vpm_a + vbm_shuffled;

	vbm_shuffled = __builtin_shuffle(vbm_shuffled, shuffle_mask_swap);
	vpm_b = __builtin_shuffle(path_metric[4+k], shuffle_mask1);
	vpm_b = vpm_b + vbm_shuffled;

	path_metric_next[2*k+1] = vpmin(vpm_a, vpm_b);

	compare1 = vpm_b < vpm_a;
	compare1 &= (v8hu){0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000};

	compare1 |= compare0;
	for (i = 0; i < 8; i++)
		path_choice |= compare1[i];

	*pc |= ((uint64_t)path_choice) << 16*k;
}

static void viterbi_feed(unsigned int start, unsigned int stop,
	uint32_t depuncture_word, uint32_t *data)
{
	unsigned int num;
	uint16_t branch_metric[8][8];
	v8hu *bm;
	v8hu *path_metric;
	v8hu *path_metric_next;
	v8hu *tmp;
	v8hu path_metric_a[8];
	v8hu path_metric_b[8];
	unsigned int jmax = 8;
	unsigned int i, j;
	uint64_t path_choice;

	num = depuncture(start, stop, depuncture_word, data, depunctured_data);

//	fprintf(stderr, "feed: %d\n", num);

	/* Special case for tail bits which are not a multiple of 32 bits */
	if (stop - start == 12) {
		num = 1;
		jmax = 6;
	}

	memcpy(path_metric_a, vb.path_metric, sizeof(path_metric_a));

	path_metric = path_metric_a;
	path_metric_next = path_metric_b;

	for (i = 0; i < num; i++) {
		viterbi_calc_branch_metric(depunctured_data[i], depuncture_word,
			branch_metric);
		for (j = 0; j < jmax; j++) {

			bm = (v8hu*)branch_metric[j];

			path_choice = 0;
			update_path_metric(bm, path_metric,
				path_metric_next, 0, &path_choice);
			update_path_metric(bm, path_metric,
				path_metric_next, 1, &path_choice);
			update_path_metric(bm, path_metric,
				path_metric_next, 2, &path_choice);
			update_path_metric(bm, path_metric,
				path_metric_next, 3, &path_choice);

			tmp = path_metric;
			path_metric = path_metric_next;
			path_metric_next = tmp;

			vb.path_choices[vb.num_paths] = path_choice;
			vb.num_paths++;
		}

#if 0
	for (int k = 0; k < 64; k++) {
		fprintf(stderr, "%d,%c", path_metric[k/8][k%8], k == 63 ? '\n' : ' ');
	}
#endif
	}

	memcpy(vb.path_metric, path_metric, sizeof(path_metric_a));
}

static __attribute__ ((noinline)) unsigned int viterbi_trace_back(uint8_t *out,
	unsigned int *len)
{
	unsigned int state = 0;
	unsigned int min_path = UINT_MAX;
	unsigned int choice;
	uint8_t out_word = 0;
	int skip = 6;
	int i;

	for (i = 0; i < 64; i++) {
		if (min_path > vb.path_metric[i]) {
			min_path = vb.path_metric[i];
			state = i;
		}
	}

	for (i = vb.num_paths-1; i >= 0; i--) {
		if (skip > 0) {
			skip--;
		} else {
			out_word |= (state & 1) << (7-(i & 7));
			if ((i & 7) == 0) {
				//printf("%d\n", out_word);
				out[i/8] = out_word;
				out_word = 0;
			}
		}
		choice = (vb.path_choices[i] >> state) & 1;
		state = (state >> 1) | (choice << 5);
	}

	if (len)
		*len = (vb.num_paths - 6) / 8;

	return min_path;
}

static void  __attribute__ ((noinline)) descramble(uint8_t *buf, unsigned int size)
{
	unsigned int state = 0x1ff;
	unsigned int feedback, next_state;
	unsigned int i, j;

	// x**9 + x**5 + 1

	for (i = 0; i < size; i++) {
		for (j = 0; j < 2; j++) {
			next_state = state << 4;
			feedback = ((state ^ next_state) >> 5) & 0xf;
			state = next_state | feedback;
		}
		buf[i] ^= state;
	}
}

static complex ofdm_symbol[2048 + 504];
static complex qpsk_symbol[1536];
static unsigned int ofdm_symbol_length;
static unsigned int frame_length;
static unsigned int frame_state;
static double phase;
static double phase_delta;

static unsigned int fft_index;
static fftw_plan fft_forward[2];
static complex fft_data[2][2048];

static uint32_t demod_data[96];

static uint32_t time_deinterleave_data[30720];
static unsigned int time_deinterleave_index;
static unsigned int time_deinterleave_ready;

static uint8_t out_data[64000];

static unsigned int subch_filter_cif_frame_start = 0;
static unsigned int subch_filter_cif_frame_end = 0;
static unsigned int subch_filter_cu_start = 0;
static unsigned int subch_filter_cu_end = 0;
static unsigned int subch_filter_num_words = 0;
static unsigned int vb_stop0, vb_stop1;
static unsigned int vb_mask0, vb_mask1;

static void write_packet(unsigned int type, unsigned int length,
	unsigned int errors, void *data)
{
	uint8_t header[8];

	if (packet_fd == -1)
		return;

	header[0] = 'L';
	header[1] = 'P';
	header[2] = 'C';
	header[3] = type;
	header[4] = (length >> 8);
	header[5] = length & 0xff;
	header[6] = (errors >> 8);
	header[7] = errors & 0xff;

	write(packet_fd, header, 8);
	write(packet_fd, data, length);
}

static unsigned int time_interleave_table[] = {
	0,
	8,
	4,
	12,
	2,
	10,
	6,
	14,
	1,
	9,
	5,
	13,
	3,
	11,
	7,
	15
};

static uint64_t frame_offset = 0;
static uint64_t frame_start = 0;
static double frame_mag = 0;
static double signal_mag = 0;

typedef double v2df __attribute__((vector_size(16)));
typedef uint64_t v2du __attribute__((vector_size(16)));
typedef double v4df __attribute__((vector_size(32)));
typedef uint64_t v4du __attribute__((vector_size(32)));

static __attribute__ ((noinline)) void demod_dqpsk(complex *a, complex *b)
{
	uint32_t *demod_data_i, *demod_data_q;
	unsigned int i, m, k;

	for (i = 1; i < 769; i++)
		qpsk_symbol[i-1] = a[i] * conj(b[i]);

	for (i = 1280; i < 2048; i++)
		qpsk_symbol[i-512] = a[i] * conj(b[i]);

	demod_data_i = demod_data;
	demod_data_q = demod_data+48;

	for (m = 0; m < 1536; ) {
		v2du mask = (v2du){0x1, 0x1};
		v2du data = (v2du){0x0, 0x0};

		for (i = 0; i < 32; i++) {
			v2df symbol;
			v2du compare;

			k = interleave_table[m];
			symbol[0] = creal(qpsk_symbol[k]);
			symbol[1] = cimag(qpsk_symbol[k]);

			compare = symbol < 0;
			data |= compare & mask;
			mask <<= 1;
			m++;
		}

		*demod_data_i++ = data[0];
		*demod_data_q++ = data[1];
	}

	if (capture_qpsk)
		write_packet(2, sizeof(qpsk_symbol), 0, qpsk_symbol);
}

static __attribute__((noinline)) void rotate(complex *data,
	float phase, float phase_delta)
{
	complex rotator;
	complex rotator_delta;
	v4df rotator_delta1;
	v4df rotator_delta2;
	v4df rotator2, r1, r2;
	v4df d, d1, d2;
	unsigned int i;
	double rotator_delta_i, rotator_delta_q;

	static const v4du ldup_mask = {0, 0, 2, 2};
	static const v4du hdup_mask = {1, 1, 3, 3};
	static const v4du swap_mask = {1, 0, 3, 2};

	rotator_delta = cosf(phase_delta) + sinf(phase_delta) * _Complex_I;
	rotator = cosf(phase) + sinf(phase) * _Complex_I;

	rotator2[0] = creal(rotator);
	rotator2[1] = cimag(rotator);
	rotator *= rotator_delta;
	rotator2[2] = creal(rotator);
	rotator2[3] = cimag(rotator);

	rotator_delta *= rotator_delta;
	rotator_delta_i = creal(rotator_delta);
	rotator_delta_q = cimag(rotator_delta);

	rotator_delta1 = (v4df){rotator_delta_i, rotator_delta_i, rotator_delta_i, rotator_delta_i};
	rotator_delta2 = (v4df){rotator_delta_q, rotator_delta_q, rotator_delta_q, rotator_delta_q};

	data += 250;

	for (i = 250; i < 2048+250; i+=2) {

		d = *(v4df *)data;
		d1 = __builtin_shuffle(d, ldup_mask); // a0 a0 a1 a1
		d2 = __builtin_shuffle(d, hdup_mask); // b0 b0 b1 b1

		d1 = d1 * rotator2;
		r1 = rotator_delta1 * rotator2;

		rotator2 = __builtin_shuffle(rotator2, swap_mask); // d0 c0 d1 c1
		d2 = d2 * rotator2;
		r2 = rotator_delta2 * rotator2;

		*(v4df *)data = _mm256_addsub_pd(d1, d2);
		rotator2 = _mm256_addsub_pd(r1, r2);

		data+=2;
	}
}

static __attribute__((noinline)) complex autocorrelate(complex *data)
{
	unsigned int i;
	complex y = 0;

	for (i = 26; i < 504-26; i++)
		y += data[i] * conj(data[i+2048]);

	return y;
}

static void decode_fic(unsigned int frame, uint32_t *data)
{
	unsigned int errors;

	switch (frame) {
	case 1:
		viterbi_reset();
		viterbi_feed(0, 2016, 0x77777777, data);
		viterbi_feed(2016, 2292, 0x37777777, data);
		viterbi_feed(2292, 2304, 0x00333333, data);
		errors = viterbi_trace_back(out_data, NULL);
		descramble(out_data, 96);
		write_packet(0, 96, errors, out_data);

		viterbi_reset();
		viterbi_feed(2304, 3072, 0x77777777, data);
		break;
	case 2:
		viterbi_feed(0, 1248, 0x77777777, data);
		viterbi_feed(1248, 1524, 0x37777777, data);
		viterbi_feed(1524, 1536, 0x00333333, data);
		errors = viterbi_trace_back(out_data, NULL);
		descramble(out_data, 96);
		write_packet(0, 96, errors, out_data);

		viterbi_reset();
		viterbi_feed(1536, 3072, 0x77777777, data);
		break;
	case 3:
		viterbi_feed(0, 480, 0x77777777, data);
		viterbi_feed(480, 756, 0x37777777, data);
		viterbi_feed(756, 768, 0x00333333, data);
		errors = viterbi_trace_back(out_data, NULL);
		descramble(out_data, 96);
		write_packet(0, 96, errors, out_data);

		viterbi_reset();
		viterbi_feed(768, 2784, 0x77777777, data);
		viterbi_feed(2784, 3060, 0x37777777, data);
		viterbi_feed(3060, 3072, 0x00333333, data);
		errors = viterbi_trace_back(out_data, NULL);
		descramble(out_data, 96);
		write_packet(0, 96, errors, out_data);
		break;
	}

//	printf("errors: %d\n", errors);
}

static void decode_msc(unsigned int cif_frame, uint32_t *data)
{
	unsigned int x = (time_deinterleave_index*subch_filter_num_words);
	unsigned int k, l;
	unsigned int y;
	uint32_t word, mask;
	uint32_t *buf;
	unsigned int cu_start = cif_frame * 48;
	unsigned int cu_end = cu_start + 48;
	unsigned int cu_size, cu_offset;
	unsigned int demod_offset = 0;
	unsigned int packet_len;
	unsigned int errors;

	if (cu_end <= subch_filter_cu_start ||
		cu_start >= subch_filter_cu_end)
		return;

	if (cu_end > subch_filter_cu_end)
		cu_end = subch_filter_cu_end;

	if (cu_start < subch_filter_cu_start) {
		demod_offset = (subch_filter_cu_start - cu_start) * 2;
		cu_start = subch_filter_cu_start;
	}

	cu_size = (cu_end - cu_start) * 8;
	cu_offset = (cu_start - subch_filter_cu_start) * 2;

	memcpy(time_deinterleave_data + x + cu_offset,
		data + demod_offset, cu_size);
	if (cu_end != subch_filter_cu_end)
		return;

	time_deinterleave_index++;
	if (time_deinterleave_index == 16) {
		time_deinterleave_index = 0;
		time_deinterleave_ready = 1;
	}

	if (!time_deinterleave_ready)
		return;

	y = time_deinterleave_index;
	for (k = 0; k < subch_filter_num_words; k++) {
		word = 0;
		mask = 0x00010001;
		for (l = 0; l < 16; l++) {
			unsigned int z = time_interleave_table[l];
			z = (y + z) % 16;
			z *= subch_filter_num_words;
			word |= time_deinterleave_data[z+k] & mask;
			mask <<= 1;
		}
		time_deinterleave_data[y*subch_filter_num_words+k] = word;
	}
	buf = &time_deinterleave_data[y*subch_filter_num_words];

	viterbi_reset();
	viterbi_feed(0, vb_stop0, vb_mask0, buf);
	viterbi_feed(vb_stop0, vb_stop1, vb_mask1, buf);
	viterbi_feed(vb_stop1, vb_stop1+12, 0x00333333, buf);
	errors = viterbi_trace_back(out_data, &packet_len);
	descramble(out_data, packet_len);
	write_packet(1, packet_len, errors, out_data);
}

static void decode_ofdm(int16_t *data, size_t size)
{
	unsigned int i;
	complex c;

	for (i = 0; i < size; i++) {
		c = data[i*2] + data[i*2+1] * _Complex_I;
		double mag = creal(c) * creal(c) + cimag(c) * cimag(c);
		signal_mag = signal_mag * 0.999999 + mag * 0.000001;

		if (frame_state <= 1) {
			frame_mag = frame_mag * 0.99 + mag * 0.01;
			if (frame_state == 0 && frame_mag < signal_mag*0.5) {
				frame_state = 1;
				frame_start = frame_offset;
			} else if (frame_state == 1 && frame_mag > signal_mag*0.6) {
				if (frame_offset - frame_start > 500) {
/*					fprintf(stderr, "ZERO: %llu %llu %llu %f!\n", frame_offset,
					frame_start, frame_offset - frame_start, delta);
					fprintf(stderr, "power: %f %f\n", signal_mag, frame_mag);*/
					frame_state = 2;
					frame_length = 0;
					phase = 0;
				} else {
					frame_state = 0;
				}
			}
		} else {
			ofdm_symbol[ofdm_symbol_length++] = c;
			if (ofdm_symbol_length == 504 + 2048) {
				unsigned int current_fft_index;
				bool do_demod = false;
				bool do_fft = false;
				unsigned int cif_frame = 0;

				if (enable_frequency_tracking) {
					complex y = 0;
					double alpha;
					y = autocorrelate(ofdm_symbol);
					alpha = carg(y) / 2048.0;
					phase_delta = phase_delta * 0.9 + alpha * 0.1;
					if (phase_delta > 3.14 / 2048.0 || phase_delta < -3.14 / 2048.0)
						phase_delta = 0;
				}

				ofdm_symbol_length = 0;

				if (frame_length < 4) {
					do_fft = true;
					if (frame_length != 0)
						do_demod = true;
				} else {
					cif_frame = (frame_length - 4) % 18;
					if (cif_frame >= subch_filter_cif_frame_start &&
					    cif_frame < subch_filter_cif_frame_end) {
						do_demod = true;
						do_fft = true;
					}

					/* The one before the first as a reference */
					if ((cif_frame+1) % 18 == subch_filter_cif_frame_start)
						do_fft = true;
				}

				if (do_fft) {
					if (enable_frequency_tracking) {
						rotate(ofdm_symbol, phase, phase_delta);
						phase += phase_delta * 2552.0;
					}
					fftw_execute(fft_forward[fft_index]);

					if (capture_fft)
						write_packet(3, sizeof(fft_data[0]), frame_length, fft_data[fft_index]);

					current_fft_index = fft_index;
					fft_index ^= 1;
				}

				if (do_demod) {
					demod_dqpsk(fft_data[current_fft_index], fft_data[fft_index]);

					if (frame_length < 4)
						decode_fic(frame_length, demod_data);
					else
						decode_msc(cif_frame, demod_data);
				}

				frame_length++;
				if (frame_length == 76)
					frame_state = 0;
			}
		}
		frame_offset++;
	}
}

void dabdecoder_init(int _read_fd, int _packet_fd)
{
	unsigned int i, k;
	int x, y;

	fft_forward[0] = fftw_plan_dft_1d(2048, ofdm_symbol+250, fft_data[0],
		FFTW_FORWARD, FFTW_ESTIMATE);
	fft_forward[1] = fftw_plan_dft_1d(2048, ofdm_symbol+250, fft_data[1],
		FFTW_FORWARD, FFTW_ESTIMATE);

	x = 0;
	k = 0;
	for (i = 0; i < 2048; i++) {
		x = (13 * x + 511) % 2048;
		if (x < 256 || x > 1792 || x == 1024)
			continue;
		y = x - 1025;
		if (y < 0)
			y += 1537;
		interleave_table[k++] = y;
	}

	read_fd = _read_fd;
	packet_fd = _packet_fd;
}

void dabdecoder_capture_qpsk(unsigned int _capture_qpsk)
{
	capture_qpsk = _capture_qpsk;
}

void dabdecoder_capture_fft(unsigned int _capture_fft)
{
	capture_fft = _capture_fft;
}

void dabdecoder_track_frequency(unsigned int _track_frequency)
{
	if (enable_frequency_tracking && !_track_frequency)
		phase_delta = 0.0;
	enable_frequency_tracking = _track_frequency;
}

int dabdecoder_process(void)
{
	static int16_t data[40000];
	int ret;

	ret = read(read_fd, data, sizeof(data));
	if (ret <= 0)
		return ret;
	decode_ofdm(data, ret / 4);

	return ret;
}

void dabdecoder_set_subch_filter(unsigned int start, unsigned int len,
	unsigned int protection_level)
{
	unsigned int n;

	subch_filter_cu_start = start;
	subch_filter_cu_end = start + len;
	subch_filter_num_words = len * 2; // 64-bits per CU

	subch_filter_cif_frame_start = start / 48;
	subch_filter_cif_frame_end = (start + len + 47) / 48;

	if (len == 0) {
		vb_stop0 = 0;
		vb_stop1 = 0;
		return;
	}

	switch (protection_level) {
	case 0: /* 1-A */
		n = len / 12;

		vb_mask0 = 0xffffffff; /* PI=24 */
		vb_mask1 = 0x7fffffff; /* PI=23 */
		vb_stop0 = 32 * (6 * n - 3);
		vb_stop1 = vb_stop0 + 31 * 3;
		break;
	case 1: /* 2-A */
		n = len / 8;

		if (n == 1) {
			vb_mask0 = 0x37373777; /* PI=13 */
			vb_mask1 = 0x37373737; /* PI=12 */
			vb_stop0 = 5 * 21;
			vb_stop1 = 20;
		} else {
			vb_mask0 = 0x37773777; /* PI=14 */
			vb_mask1 = 0x37373777; /* PI=13 */
			vb_stop0 = 22 * (2*n - 3);
			vb_stop1 = vb_stop0 + 21 * (4*n + 3);
		}

		break;
	case 2: /* 3-A */
		n = len / 6;

		vb_mask0 = 0x33333333; /* PI=8 */
		vb_mask1 = 0x13333333; /* PI=7 */
		vb_stop0 = 16 * (6 * n - 3);
		vb_stop1 = vb_stop0 + 15 * 3;
		break;
	case 3: /* 4-A */
		n = len / 4;

		vb_mask0 = 0x11131313; /* PI=3 */
		vb_mask1 = 0x11131113; /* PI=2 */
		vb_stop0 = 11 * (4 * n - 3);
		vb_stop1 = vb_stop0 + 10 * (2 * n + 3);
		break;
	case 4: /* 1-B */
		n = len / 27;

		vb_mask0 = 0x33373337; /* PI=10 */
		vb_mask1 = 0x33333337; /* PI=9 */
		vb_stop0 = 18 * (24 * n - 3);
		vb_stop1 = vb_stop0 + 17 * 3;
		break;
	case 5: /* 2-B */
		n = len / 21;

		vb_mask0 = 0x13331333; /* PI=6 */
		vb_mask1 = 0x13131333; /* PI=5 */
		vb_stop0 = 14 * (24 * n - 3);
		vb_stop1 = vb_stop0 + 13 * 3;
		break;
	case 6: /* 3-B */
		n = len / 18;

		vb_mask0 = 0x13131313; /* PI=4 */
		vb_mask1 = 0x11131313; /* PI=3 */
		vb_stop0 = 12 * (24 * n - 3);
		vb_stop1 = vb_stop0 + 11 * 3;
		break;
	case 7: /* 4-B */
		n = len / 15;

		vb_mask0 = 0x11131113; /* PI=2 */
		vb_mask0 = 0x11111113; /* PI=1 */
		vb_stop0 = 10 * (24 * n - 3);
		vb_stop1 = vb_stop0 + 9 * 3;
		break;
	}

	vb_stop0 *= 4;
	vb_stop1 *= 4;

	time_deinterleave_index = 0;
	time_deinterleave_ready = 0;
}

int main(void)
{
	dabdecoder_init(0, 1);
	dabdecoder_set_subch_filter(0, 60, 2);

	while (true)
		if (dabdecoder_process() <= 0)
			break;

	return 0;
}
