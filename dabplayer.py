import gi
from gi.repository import Gst, GObject
import os

class DABAudioPlayer:

	def __init__(self):
		self.src_idle_id = None
		self.frames = []
		self.do_need_data = False
		self.pipeline = None
		self.standby = False

	def create_pipeline(self):
		self.pipeline = Gst.Pipeline.new('dabplayer')
		self.src = Gst.ElementFactory.make('appsrc', 'source')
		aacparse = Gst.ElementFactory.make('aacparse', 'aacparse')
		avdec_aac = Gst.ElementFactory.make('avdec_aac', 'avdec_aac')
		audioconvert = Gst.ElementFactory.make('audioconvert', 'audioconvert')
		audioresample = Gst.ElementFactory.make('audioresample', 'audioresample')
		sink = Gst.ElementFactory.make('pulsesink', 'pulseink')

		self.pipeline.add(self.src)
		self.pipeline.add(aacparse)
		self.pipeline.add(avdec_aac)
		self.pipeline.add(audioconvert)
		self.pipeline.add(audioresample)
		self.pipeline.add(sink)

		self.src.link(aacparse)
		aacparse.link(avdec_aac)
		avdec_aac.link(audioconvert)
		audioconvert.link(audioresample)
		audioresample.link(sink)

		self.src.connect('need-data', self.need_data)
		self.src.connect('enough-data', self.enough_data)

	def push_frame(self, frame):
		self.frames.append(frame)
		if len(self.frames) > 4:
			if self.standby:
				self.pipeline.set_state(Gst.State.PLAYING)
				self.standby = False
			if self.do_need_data and self.src_idle_id is None:
				self.src_idle_id = GObject.idle_add(self.feed_data)

	def start_playback(self):
		if self.pipeline:
			self.stop_playback()
		self.create_pipeline()
		self.standby = True

	def stop_playback(self):
		if self.pipeline is None:
			return
		self.pipeline.set_state(Gst.State.NULL)
		if self.src_idle_id is not None:
			GObject.source_remove(self.src_idle_id)
			self.src_idle_id = None
		self.frames = []
		self.pipeline = None

	def feed_data(self):
		try:
			frame = self.frames.pop(0)
		except:
			self.src_idle_id = None
			return False

		buf = Gst.Buffer.new_wrapped(frame)
		self.src.emit('push-buffer', buf)

		return True

	def need_data(self, src, size):
		if self.src_idle_id is None:
			self.src_idle_id = GObject.idle_add(self.feed_data)
			self.do_need_data = True

	def enough_data(self, src):
		self.do_need_data = False
		if self.src_idle_id is not None:
			GObject.source_remove(self.src_idle_id)
			self.src_idle_id = None

Gst.init(None)
