#!/usr/bin/env python3

import struct
import crc16
import sys
import os
import dabplayer
import threading
import time

class DABDatabase:

	def __init__(self):
		self.ensemble = None
		self.service_organizations = {}
		self.sub_channels = {}
		self.labels = {}
		self.programme_type = {}
		self.global_defs = {}
		self.packet_mode = {}

	def __update_dict(self, d, key, obj):
		old = d.get(key)
		if old != obj:
			print('Update', old, obj)
			d[key] = obj
			return True
		return False

	def update(self, obj):
		changed = False
		t = type(obj)

		if t is DABEnsemble:
			if self.ensemble:
				if self.ensemble.eid != obj.eid:
					self.labels = {}
					self.sub_channels = {}
					self.service_organizations = {}
					self.programme_type = {}
					self.global_defs = {}
					self.packet_mode = {}
					changed = True
			else:
				changed = True
			self.ensemble = obj
		elif t is DABLabel:
			changed = self.__update_dict(self.labels, (obj.ext, obj.lid), obj)
		elif t is DABSubChannel:
			changed = self.__update_dict(self.sub_channels, obj.subchid, obj)
		elif t is DABServiceOrganization:
			changed = self.__update_dict(self.service_organizations, obj.sid, obj)
		elif t is DABProgrammeType:
			changed = self.__update_dict(self.programme_type, obj.sid, obj)
		elif t is DABServiceComponentGlobalDefinition:
			changed = self.__update_dict(self.global_defs, (obj.sid, obj.subchid, obj.scid), obj)
		elif t is DABServiceComponentPacketMode:
			changed = self.__update_dict(self.packet_mode, obj.scid, obj)

		return changed

	def __get_secondary_label_index(self, gd_index):
		gd = self.global_defs.get(gd_index)
		if gd is None:
			return None
		return (4, (gd.sid, gd.scids))

	def get_label(self, obj):
		t = type(obj)
		if t is DABEnsemble:
			index = (0, obj.eid)
		elif t is DABServiceOrganization:
			index = (1, obj.sid)
		elif t is DABServiceComponentAudioStream:
			if obj.primary:
				index = (1, obj.sid)
			else:
				index = self.__get_secondary_label_index((obj.sid, obj.subchid, None))
				if index is None:
					return None
		elif t is DABServiceComponentPacketData:
			if obj.primary:
				index = (5, obj.sid)
			else:
				index = self.__get_secondary_label_index((obj.sid, None, obj.scid))
				if index is None:
					return None
		else:
			return None

		if index in self.labels:
			return self.labels[index].text.strip()
		else:
			return None

	def get_sub_channel(self, sc):
		t = type(sc)
		if t is DABServiceComponentAudioStream:
			return self.sub_channels.get(sc.subchid)
		elif t is DABServiceComponentPacketData:
			pm = self.packet_mode[sc.scid]
			return self.sub_channels.get(pm.subchid)
		else:
			return None

class DABError:
	def __repr__(self):
		return 'DABError'

class DABErrorFIBIncorrectCRC(DABError):
	def __repr__(self):
		return 'DABErrorFIBIncorrectCRC'

class DABErrorAUIncorrectCRC(DABError):
	pass

class DABObject:
	def __eq__(self, other):
		if isinstance(other, self.__class__):
			return self.__dict__ == other.__dict__
		return NotImplemented

	def __ne__(self, other):
		if isinstance(other, self.__class__):
			return self.__dict__ != other.__dict__
		return NotImplemented

class DABEnsemble(DABObject):

	def __init__(self, eid, change_flag, alarm, cif_count, occurrence_change):
		self.eid = eid
		self.change_flag = change_flag
		self.alarm = alarm
		self.cif_count = cif_count
		self.occurrence_change = occurrence_change

	def __repr__(self):
		return 'DABEnsemble(id={}, change_flag={}, alarm={}, cif_count={}, occurrence_change={})'.format(self.eid, self.change_flag, self.alarm, self.cif_count, self.occurrence_change)

class DABLabel(DABObject):

	def __init__(self, ext, lid, charset, text):
		self.ext = ext
		self.lid = lid
		self.charset = charset
		self.text = text

	def __repr__(self):
		return 'DABLabel(ext={}, id={}, charset={}, text=\'{}\')'.format(self.ext, self.lid, self.charset, self.text)

class DABSubChannel(DABObject):

	def __init__(self, subchid, startaddr, size, option, protection_level):
		self.subchid = subchid
		self.startaddr = startaddr
		self.size = size
		self.option = option
		self.protection_level = protection_level

	def __repr__(self):
		return 'DABSubChannel(id={}, startaddr={}, size={}, option={}, protection_level={})'.format(self.subchid, self.startaddr, self.size, self.option, self.protection_level)

class DABServiceComponent(DABObject):

	def __init__(self, tmid, primary, ca):
		self.tmid = tmid
		self.primary = primary
		self.ca = ca

	def __repr__(self):
		return 'DABServiceComponent(tmid={}, primary={}, ca={})'.format(self.tmid, self.primary, self.ca)

class DABServiceComponentAudioStream(DABServiceComponent):

	def __init__(self, tmid, primary, ca, ascty, subchid):
		super(DABServiceComponentAudioStream, self).__init__(tmid, primary, ca)
		self.ascty = ascty
		self.subchid = subchid

	def __repr__(self):
		return 'DABServiceComponentAudioStream(tmid={}, ascty={}, subchid={}, primary={}, ca={})'.format(self.tmid, self.ascty, self.subchid, self.primary, self.ca)

class DABServiceComponentPacketData(DABServiceComponent):

	def __init__(self, tmid, primary, ca, scid):
		super(DABServiceComponentPacketData, self).__init__(tmid, primary, ca)
		self.scid = scid

	def __repr__(self):
		return 'DABServiceComponentPacketData(tmid={}, scid={}, primary={}, ca={})'.format(self.tmid, self.scid, self.primary, self.ca)

class DABServiceOrganization(DABObject):

	def __init__(self, sid, local, caid):
		self.sid = sid
		self.local = local
		self.caid = caid
		self.components = []

	def __repr__(self):
		return 'DABServiceOrganization(id={}, local={}, caid={})'.format(self.sid, self.local, self.caid)

	def add_component(self, c):
		self.components.append(c)

class DABDateTime(DABObject):

	def __init__(self, lsi, mjd, hours, minutes, seconds, milliseconds):
		self.lsi = lsi
		self.mjd = mjd
		self.hours = hours
		self.minutes = minutes
		self.seconds = seconds
		self.milliseconds = milliseconds

	def __repr__(self):
		return 'DABDateTime(lsi={}, date={}, time={:02d}:{:02d}:{:02d}:{:03d}'.format(self.lsi, self.mjd, self.hours, self.minutes, self.seconds, self.milliseconds)

class DABProgrammeType(DABObject):

	def __init__(self, sid, dynamic, code):
		self.sid = sid
		self.dynamic = dynamic
		self.code = code

	def __repr__(self):
		return 'DABProgrammeType(sid={}, dynamic={}, code={}'.format(self.sid, self.dynamic, self.code)

class DABCountryLTOInternationalTable(DABObject):

	def __init__(self, lto, ecc, international_table):
		self.lto = lto
		self.ecc = ecc
		self.international_table = international_table

	def __repr__(self):
		return 'DABCountryLTOInternationalTable(lto={}, ecc={}, internation_table={}'.format(self.lto, self.ecc, self.intrnation_table)

class DABServiceComponentGlobalDefinition(DABObject):

	def __init__(self, sid, scids, scid, subchid):
		self.sid = sid
		self.scids = scids
		self.scid = scid
		self.subchid = subchid

	def __repr__(self):
		return 'DABServiceComponentGlobalDefinition(sid={}, scids={}, subchid={}, scid={})'.format(self.sid, self.scids, self.subchid, self.scid)

class DABServiceComponentPacketMode(DABObject):

	def __init__(self, scid, use_data_groups, dscty, subchid, packet_address, caorg):
		self.scid = scid
		self.use_data_groups = use_data_groups
		self.dscty = dscty
		self.subchid = subchid
		self.packet_address = packet_address
		self.caorg = caorg

	def __repr__(self):
		return 'DABServiceComponentPacketMode(scid={}, use_data_groups={}, dscty={}, subchid={}, packet_address={}, caorg={})'.format(self.scid , self.use_data_groups, self.dscty, self.subchid, self.packet_address, self.caorg)

class DABFICParser:

	def __init__(self):
		pass

	def read_int(self, data, offset, l):
		val = 0
		for i in range(0, l):
			val <<= 8
			val |= data[offset+i]
		return val

	def handle_fig0_ext0(self, fig, fig_len):
		eid = self.read_int(fig, 1, 2)
		cf = fig[3] >> 6
		alarm = (fig[3] >> 5) & 1
		cif_count = self.read_int(fig, 1, 2) & 0x1fff
#		Some stations don't send OC?
#		oc = fig[6]
		oc = 0
		yield DABEnsemble(eid, cf, alarm, cif_count, oc)

	def handle_fig0_ext1(self, fig, fig_len):
		x = 1
		while x < fig_len:
			subchid = fig[x] >> 2
			startaddr = self.read_int(fig, x, 2) & 0x3ff
			lng = fig[x+2] >> 7
			opt = (fig[x+2] >> 4) & 0x7
			pl = (fig[x+2] >> 2) & 0x3
			size = self.read_int(fig, x+2, 2) & 0x3ff
			yield DABSubChannel(subchid, startaddr, size, opt, pl)
			x += 3 + lng

	def handle_fig0_ext2(self, fig, fig_len, pd):
		x = 1
		while x < fig_len:
			if pd:
				id_len = 4
			else:
				id_len = 2
			sid = self.read_int(fig, x, id_len)
			x += id_len
			local = fig[x] >> 7
			caid = (fig[x] >> 4) & 0x7
			nsc = fig[x] & 0xf
			so = DABServiceOrganization(sid, local, caid)
			x += 1
			for i in range(0, nsc):
				tmid = fig[x] >> 6
				primary = (fig[x+1] >> 1) & 1
				ca = fig[x+1] & 1

				if tmid == 0:
					ascty = fig[x] & 0x3f
					subchid = fig[x+1] >> 2
					sc = DABServiceComponentAudioStream(tmid, primary, ca, ascty, subchid)
				elif tmid == 3:
					scid = (self.read_int(fig, x, 2) >> 2) & 0x3ff
					sc = DABServiceComponentPacketData(tmid, primary, ca, scid)
				else:
					sc = DABServiceComponent(tmid, primary, ca)
				sc.sid = sid
				so.add_component(sc)
				x += 2
			yield so

	def handle_fig0_ext3(self, fig, fig_len):
		x = 1
		while x < fig_len:
			scid = self.read_int(fig, x, 2) >> 4
			caorg_flag = fig[x+1] & 1
			dg = fig[x+2] >> 7
			dscty = fig[x+2] & 0x3f
			subchid = fig[x+3] >> 2
			packet_address = self.read_int(fig, x+3, 2) & 0x3ff
			x += 5
			if caorg_flag:
				caorg = self.read_int(fig, x, 2)
				x += 2
			else:
				caorg = None
			yield DABServiceComponentPacketMode(scid, dg, dscty, subchid, packet_address, caorg)

	def handle_fig0_ext5(self, fig, fig_len):
		x = 1
		while x < fig_len:
			ls = fig[x] >> 7
			if ls:
				scid = self.read_int(x, 2) & 0x3ff
				subchid = None
				x += 2
			else:
				subchid = fig[x] & 0x3f
				scid = None
				x += 1
			language = fig[x]
			x += 1

	def handle_fig0_ext8(self, fig, fig_len, pd):
		x = 1
		while x < fig_len:
			if pd:
				sid = self.read_int(fig, x, 4)
				x += 4
			else:
				sid = self.read_int(fig, x, 2)
				x += 2
			ext = fig[x] >> 7
			scids = fig[x] & 0xf

			ls = fig[x+1] >> 7
			if ls:
				scid = self.read_int(fig, x+1, 2) & 0x3ff
				subchid = None
			else:
				subchid = fig[x+1] & 0x3f
				scid = None

			yield DABServiceComponentGlobalDefinition(sid, scids, scid, subchid)

			x += 2 + ls + ext

	def handle_fig0_ext9(self, fig, fig_len):
		ext = fig[1] >> 7
		lto = fig[1] & 0x3f
		ecc = fig[2]
		table = fig[3]

		yield DABCountryLTOInternationalTable(lto, ecc, table)

	def handle_fig0_ext10(self, fig, fig_len):
		mjd = (self.read_int(fig, 1, 3) >> 6) & 0x1ffff
		lsi = (fig[3] >> 5) & 1
		utc = (fig[3] >> 3) & 1
		time = self.read_int(fig, 3, 2)
		hours = (time >> 6) & 0x1f
		minutes = time & 0x3f
		if utc:
			time = self.read_int(fig, 5, 2)
			seconds = (time >> 10) & 0x3f
			milliseconds = time & 0x3ff
		else:
			seconds = 0
			milliseconds = 0

		yield DABDateTime(lsi, mjd, hours, minutes, seconds, milliseconds)


	def handle_fig0_ext17(self, fig, fig_len):
		x = 1
		while x < fig_len:
			sid = self.read_int(fig, x, 2)
			dynamic = fig[x+2] >> 7
			code = fig[x+3] & 0x1f
			yield DABProgrammeType(sid, dynamic, code)
			x += 4

	def handle_fig(self, t, l, fig):
		if t == 0x1:
			charset = fig[0] >> 4
			ext = fig[0] & 0x7

			if ext == 4 or ext == 6:
				pd = fig[1] >> 7
				scsid = fig[1] & 0xf
				if pd:
					lid = self.read_int(fig, 2, 4)
					x = 6
				else:
					lid = self.read_int(fig, 2, 2)
					x = 4
				if ext == 6:
					app_type = fig[x] & 0x1f
					lid = (lid, scsid, app_type)
					x += 1
				else:
					lid = (lid, scsid)
			elif ext == 0 or ext == 1:
				lid = self.read_int(fig, 1, 2)
				x = 3
			elif ext == 5:
				lid = self.read_int(fig, 1, 4)
				x = 5
			else:
				print('Unknown Label type')
				return

			text = ''.join(map(chr, fig[x:x+16]))
			yield DABLabel(ext, lid, charset, text)
		elif t == 0x0:
			cn = fig[0] >> 7
			oe = (fig[0] >> 6) & 1
			pd = (fig[0] >> 5) & 1
			ext = fig[0] & 0x1f
			if ext == 0:
				yield from self.handle_fig0_ext0(fig, l)
			elif ext == 1:
				yield from self.handle_fig0_ext1(fig, l)
			elif ext == 2:
				yield from self.handle_fig0_ext2(fig, l, pd)
			elif ext == 3:
				yield from self.handle_fig0_ext3(fig, l)
			elif ext == 5:
#				yield from self.handle_fig0_ext5(fig, l)
				pass
			elif ext == 8:
				yield from self.handle_fig0_ext8(fig, l, pd)
			elif ext == 9:
				yield from self.handle_fig0_ext9(fig, l)
			elif ext == 10:
				yield from self.handle_fig0_ext10(fig, l)
			elif ext == 17:
				yield from self.handle_fig0_ext17(fig, l)
			else:
#				print('MCI and SI:', ext, cn, oe, pd, l)
				pass
		else:
			print('Unknown FIG', t, l)
			pass

	def handle_fib(self, fib):
		crc = crc16.crc16xmodem(fib, 0xffff)
		if crc != 0x1d0f:
			return DABErrorFIBIncorrectCRC()
		x = 0
		while x < 30:
			if fib[x] == 0xff:
				break
			l = fib[x] & 0x1f
			t = fib[x] >> 5

			if x + l > 30:
				print('Warning: Too long FIB', l, x+l)
				break

#			try:
			yield from self.handle_fig(t, l, fib[x+1:x+l+1])
#			except:
#				yield DABError()
			x += 1 + l

	def handle_fic(self, data):
		for i in range(0, len(data), 32):
			yield from self.handle_fib(data[i:i+32])

class DABPlusAudioParser:

	def __init__(self):
		self.frames = []
		self.sync = False

	def generate_audio_frames(self):

		dac_rate = (self.superframe[2] >> 6) & 1
		sbr_flag = (self.superframe[2] >> 5) & 1
		num_aus = 2 + dac_rate
		if sbr_flag == 0:
			num_aus *= 2

		j = 3
		au_start = [3 + ((num_aus-1)*12+4)//8]
		for i in range(1, num_aus):
			x = (self.superframe[j] << 8) | self.superframe[j+1]
			if i % 2 == 1:
				x >>= 4
				j += 1
			else:
				j += 2
			x &= 0xfff
			au_start.append(x)
		au_start.append(len(self.superframe)//120*110)


		if dac_rate:
			if sbr_flag:
				adts_rate = 6
			else:
				adts_rate = 3
		else:
			if sbr_flag:
				adts_rate = 8
			else:
				adts_rate = 5

		for i in range(0, len(au_start)-1):
			au_frame = self.superframe[au_start[i]:au_start[i+1]]
			crc = crc16.crc16xmodem(au_frame, 0xffff)
			if crc != 0x1d0f:
				print('AU: Incorrect CRC {:x}'.format(crc))
				continue

			au_size = au_start[i+1] - au_start[i] - 2
			adts_size = au_size + 7
			header = [0] * 7
			header[0] = 0xFF
			header[1] = 0xF0
			header[1] |= 0x08
			header[1] |= 0x01
			header[2] |= adts_rate<<2;
			header[3] |= 0x08
			header[3] |= 0x04
			header[3] |= 2<<6
			header[3] |= (adts_size >> 11) & 0xff;
			header[4] = (adts_size >> 3) & 0xff;
			header[5] |= (adts_size & 0x07) << 5;
			header[5] |= 0x1F
			header[6] |= 0xFC
			frame = bytes(header)
			frame += au_frame[:-2] # Skip CRC
			yield frame


	def firecode(self, data):
		crc = 0x0000

		for d in data:
			d = d << 9
			for i in range(0, 8):
				crc <<= 1
				if (d ^ crc) & 0x10000:
					crc ^= 0x782f
				d <<= 1
		return crc & 0xffff

	def check_superframe(self):
		crc = self.firecode(self.superframe[2:11]+self.superframe[0:2])
		return crc == 0x0000

	def handle_cif(self, data):
		self.frames.append(data)
		if len(self.frames) > 5:
			self.frames.pop(0)
		elif len(self.frames) < 5:
			return

		self.superframe = bytes()
		for f in self.frames:
			self.superframe += f

		if not self.check_superframe():
			if self.sync:
				print('Desync')
				self.sync = False
			return

		self.sync = True

		yield from self.generate_audio_frames()
		self.frames = []

	def reset(self):
		self.sync = False
		self.frames = []
