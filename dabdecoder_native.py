#!/usr/bin/env python3

from ctypes import Structure, c_char_p, c_uint, c_int, c_size_t, \
		c_ssize_t, c_char, c_void_p, c_bool, create_string_buffer, \
		POINTER as _POINTER, CDLL as _cdll, memmove as _memmove, byref as _byref

_lib = _cdll('./dabdecoder.so')

init = _lib.dabdecoder_init
init.argtypes = (c_int, c_int)

process = _lib.dabdecoder_process
process.restypes = c_int

capture_qpsk = _lib.dabdecoder_capture_qpsk
capture_qpsk.argtypes = (c_int, )

capture_fft = _lib.dabdecoder_capture_fft
capture_fft.argtypes = (c_int, )

track_frequency = _lib.dabdecoder_track_frequency
track_frequency.argtypes = (c_int, )

set_subch_filter = _lib.dabdecoder_set_subch_filter
set_subch_filter.argtypes = (c_int, c_int, c_int)
