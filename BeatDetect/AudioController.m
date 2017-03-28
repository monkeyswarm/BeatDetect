#import "AudioController.h"

#import <AVFoundation/AVFoundation.h>
#import <AudioToolbox/AudioQueue.h>

#import "BDPAudioUtils.h"

#define OVERLAP_FRAMECOUNT 200 // Duration, in samples, of the crossfade between playback chunks.
#define NUM_CHANNELS 2 // Decode to stereo.
#define NUM_BUFFERS 3 // Use three buffers for audio queue playback.
#define BUFFER_SIZE 4096
#define SAMPLE_RATE 44100
#define CLICK_FRAMECOUNT 1000 // Duration, in samples, of the click sound.
/// Multiplier for the click audio. This should be in terms of SAMPLE_TYPE, i.e. a little less than
// 2^15 for int16_t, a little less than +1.0 for float.
#define CLICK_AMPLITUDE 32000

@interface AudioController()
@end


@implementation AudioController {
  // Ramp which goes from 0.0 to 1.0 every beat.
  CGFloat _phasor;
  // Time duration of one beat. If zero, then play audio normal speed
  NSTimeInterval _beatInterval;
  // Index into the audio data.
  NSUInteger _frameIndex;

  int16_t* _audioData;
  NSUInteger _audioDataFrameCount;

  // Array of long, representing sample offset for each subdivision. This means it has
  // 4 entries for each beat.
  long *_beatMap;
  int _beatMapCount;

  // The currently playing beat index.
  NSUInteger _currBeat;
  // The current subdivision, 0 to 3, within the current beat.
  NSUInteger _currSubdiv;

  // Bookkeeping for playback.
  NSUInteger _trackSampleStartInPrevSubdiv;
  NSUInteger _trackSampleOffsetInPrevSubdiv;
  NSUInteger _trackSampleStartInSubdiv;
  NSUInteger _trackSampleOffsetInSubdiv;

  // x-fade amplitude table, for xfading between each subdivision to avoid clicks. Equal power xfade
  // so max val is .7071.
  CGFloat _overlapTable[OVERLAP_FRAMECOUNT];
  // Envelope table (ramp from 1 to 0) for clicks.
  CGFloat _clickEnvelopeTable[CLICK_FRAMECOUNT];

  // Bookkeeping for clicking.
  NSUInteger _clickStart; // Sample offset index of current click.
  BOOL _isInClick; // Whether we are currently rendering a click.

  AudioQueueRef _queue;
}

void callback(void *custom_data, AudioQueueRef queue, AudioQueueBufferRef buffer) {
  AudioController *controller = (__bridge AudioController *)custom_data;
  if (!controller) {
    return;
  }

  // If _audioData has not yet been allocated (as when this is called in setupAQ), pass back empty
  // buffers.
  if (!controller->_audioData) {
    memset(buffer->mAudioData, 0, buffer->mAudioDataByteSize); // Clear buffer.
    AudioQueueEnqueueBuffer(queue, buffer, 0, NULL);
    return;
  }

  int16_t* casted_buffer = (int16_t*)buffer->mAudioData;
  memset(casted_buffer, 0, buffer->mAudioDataByteSize); // Clear buffer.

  NSUInteger inNumberFrames = buffer->mAudioDataByteSize / sizeof(int16_t) / NUM_CHANNELS;
  BOOL finished = NO;
  if (controller->_beatInterval == 0) { // Normal playback.
    for (NSUInteger i = 0; i < inNumberFrames * NUM_CHANNELS; i += NUM_CHANNELS ) {
      for (int j=0; j<NUM_CHANNELS; j++) {
        casted_buffer[i+j] = controller->_audioData[controller->_frameIndex * NUM_CHANNELS + j];
      }

      // Check sample position against beat offsets to send delegate call and render click.
      if (controller->_currBeat * 4 < controller->_beatMapCount) {
        NSUInteger currBeatFrameOffset = controller->_beatMap[controller->_currBeat * 4] ;
        // Check if playback has progressed to the next beat.
        if (controller->_frameIndex > currBeatFrameOffset) {
          controller->_isInClick = YES;
          controller->_clickStart = controller->_frameIndex;
          controller->_currBeat++;
          [controller beatDidChange];
        }
        // Render click.
        if (controller->_playbackClicks && controller->_isInClick) {
          NSUInteger framesSinceClickStart = controller->_frameIndex - controller->_clickStart;
          if (framesSinceClickStart < CLICK_FRAMECOUNT) {
            // Compute amplitude, apply to noise.
            CGFloat amp = controller->_clickEnvelopeTable[framesSinceClickStart];
            casted_buffer[i] = casted_buffer[i+1] = drand48() * amp * CLICK_AMPLITUDE;
          } else {
            controller->_isInClick = NO;
          }
        }
      }
      // Increment the current frame in the controller.
      controller->_frameIndex++;
    }
  } else { // BPM-based playback.
    // The amount that _phasor will be increased on each frame.
    CGFloat phasorIncrement = (1.0 / controller->_beatInterval) / SAMPLE_RATE;

    for (NSUInteger i = 0; i < inNumberFrames * NUM_CHANNELS; i+=NUM_CHANNELS ) {
      // get quarter of a beat, using phasor
      NSUInteger subdiv = (NSUInteger)(controller->_phasor * 4); // 0-3
      if (controller->_currSubdiv != subdiv) { // New subdivision of a beat.
        controller->_currSubdiv = subdiv;
        if (subdiv == 0) { // New beat.
          controller->_isInClick = YES;
          controller->_clickStart = controller->_frameIndex;
          controller->_currBeat++;
          [controller beatDidChange];
        }

        NSUInteger beatMapIndex = controller->_currBeat * 4 + controller->_currSubdiv;
        // check if still can use beat map
        if (beatMapIndex >= controller->_beatMapCount) {
          finished = YES;
          break;
        }

        // Move current subdivision bookmark values to previous subdivision bookmarks.
        controller->_trackSampleStartInPrevSubdiv = controller->_trackSampleStartInSubdiv;
        controller->_trackSampleOffsetInPrevSubdiv = controller->_trackSampleOffsetInSubdiv;
        // Get new bookmark values for the subdivision we will start playing.
        controller->_trackSampleOffsetInSubdiv = 0;
        controller->_trackSampleStartInSubdiv = controller->_beatMap[beatMapIndex];
      }

      // Copy and/or xfade samples from audio data into the output buffer.

      // Within the overlap/xfade between subdivisions.
      if(controller->_trackSampleOffsetInSubdiv < OVERLAP_FRAMECOUNT) {
        for (NSUInteger j = 0; j < NUM_CHANNELS; j++) {
          NSUInteger sampleCurrIndex =
              controller->_trackSampleStartInSubdiv + controller->_trackSampleOffsetInSubdiv;
          int16_t sampleCurrVal = controller->_audioData[sampleCurrIndex * NUM_CHANNELS + j];

          NSUInteger samplePrevIndex =
              controller->_trackSampleStartInPrevSubdiv + controller->_trackSampleOffsetInPrevSubdiv;
          int16_t samplePrevVal = controller->_audioData[samplePrevIndex * NUM_CHANNELS + j];
          // xfade previous and current sample values.
          CGFloat xfadeAmpCurr = controller->_overlapTable[controller->_trackSampleOffsetInSubdiv];
          CGFloat xfadeAmpPrev =
              controller->_overlapTable[OVERLAP_FRAMECOUNT - 1 - controller->_trackSampleOffsetInSubdiv];
          int16_t sampleOutputVal = sampleCurrVal * xfadeAmpCurr + samplePrevVal * xfadeAmpPrev;
          casted_buffer[i+j] = sampleOutputVal;
        }
        controller->_trackSampleOffsetInSubdiv++;
        controller->_trackSampleOffsetInPrevSubdiv++;
      } else { // Not within the xfade between subdivisions.
        for (NSUInteger j = 0; j < NUM_CHANNELS; j++) {
          NSUInteger sampleIndex =
              controller->_trackSampleStartInSubdiv + controller->_trackSampleOffsetInSubdiv;
          int16_t sampleOutputVal =
              controller->_audioData[sampleIndex * NUM_CHANNELS + j] * .7071; // Same max val as the equal-power xfade.
          casted_buffer[i+j] = sampleOutputVal;
        }
        controller->_trackSampleOffsetInSubdiv++;
      }

      // Render click. DEI combine with previous copy.
      if (controller->_playbackClicks && controller->_isInClick) {
        NSUInteger framesSinceClickStart = controller->_frameIndex - controller->_clickStart;
        if (framesSinceClickStart < CLICK_FRAMECOUNT) {
          // Compute amplitude, apply to noise.
          CGFloat amp = controller->_clickEnvelopeTable[framesSinceClickStart];
          casted_buffer[i] = casted_buffer[i+1] = drand48() * amp * CLICK_AMPLITUDE;
        } else {
          controller->_isInClick = NO;
        }
      }

      // Increment / wrap phasor 0-1.
      controller->_phasor += phasorIncrement;
      if (controller->_phasor >= 1.0) {
        controller->_phasor -= 1.0;
      }

      controller->_frameIndex++;
    }
  }

  AudioQueueEnqueueBuffer(queue, buffer, 0, NULL);

  // See if we're done.
  if (controller->_frameIndex > controller->_audioDataFrameCount || finished) {
    [controller playbackDidFinish];
  }
}

- (instancetype)init {
  self = [super init];
  if (self) {
    // Make overlap table.
    for(NSUInteger i = 0; i < OVERLAP_FRAMECOUNT; i++){
      _overlapTable[i] = (CGFloat)sqrt((CGFloat)i/OVERLAP_FRAMECOUNT) * .7071;
    }
    // Make click ramp table.
    for (NSUInteger i = 0; i < CLICK_FRAMECOUNT; i++) {
      _clickEnvelopeTable[i] = (CLICK_FRAMECOUNT - i) / (CGFloat)CLICK_FRAMECOUNT;
    }

    [self setupAQ];
  }
  return self;
}

- (void)dealloc{
  AudioQueueStop(_queue, false);
  AudioQueueDispose(_queue, false);
}

- (void)setBPM:(CGFloat)BPM {
  if (BPM == 0) {
    _beatInterval = 0;
  } else {
    _beatInterval = 60 / BPM;
  }
}

- (void)jumpToBeatIndex:(NSUInteger)beatIndex {
  _currBeat = beatIndex;
  if (_beatMapCount > _currBeat * 4) {
    _frameIndex = _beatMap[_currBeat * 4];
  } else {
    // DEI error, out of bounds
    _frameIndex = 0;
  }
}

- (void)stop {
  AudioQueueStop(_queue, false);
  AudioQueueDispose(_queue, false);

  if (_audioData) {
    free(_audioData);
    _audioData = nil;
  }
  if (_beatMap) {
    free(_beatMap);
    _beatMap = nil;
  }
  _beatMapCount = 0;
  _frameIndex = 0;
  _audioDataFrameCount = 0;;
  _currBeat = 0;
  _currSubdiv = 0;
  _trackSampleStartInPrevSubdiv = 0;
  _trackSampleOffsetInPrevSubdiv = 0;
  _trackSampleStartInSubdiv = 0;
  _trackSampleOffsetInSubdiv = 0;
}

- (void)setupAQ {
  AudioStreamBasicDescription format;
  format.mSampleRate       = SAMPLE_RATE;
  format.mFormatID         = kAudioFormatLinearPCM;
  format.mFormatFlags      = kLinearPCMFormatFlagIsSignedInteger | kAudioFormatFlagIsPacked;
  format.mBitsPerChannel   = 8 * sizeof(int16_t);
  format.mChannelsPerFrame = NUM_CHANNELS;
  format.mBytesPerFrame    = sizeof(int16_t) * NUM_CHANNELS;
  format.mFramesPerPacket  = 1;
  format.mBytesPerPacket   = format.mBytesPerFrame * format.mFramesPerPacket;
  format.mReserved         = 0;
  AudioQueueBufferRef buffers[NUM_BUFFERS];

  void *context = (__bridge void *)self;
  AudioQueueNewOutput(&format, callback, context, /*CFRunLoopGetCurrent()*/nil, nil, 0, &_queue);

  for (unsigned int i = 0; i < NUM_BUFFERS; i++) {
    AudioQueueAllocateBuffer(_queue, BUFFER_SIZE, &buffers[i]);
    buffers[i]->mAudioDataByteSize = BUFFER_SIZE;
    callback(context, _queue, buffers[i]);
  }
}

- (void)play {
  AudioQueueStart(_queue, NULL);
}

// Chop 1-entry-per-beat into a list of 4 subdivisions per beat.
- (NSArray *)chopBeatMap:(NSArray *)beatMap {
  if (beatMap.count < 1) {
    return nil;
  }

  NSMutableArray *choppedBeatMap = [NSMutableArray array];
  for (NSUInteger beatIndex = 0; beatIndex < ([beatMap count] - 1); beatIndex++) {
    double beatA = [beatMap[beatIndex] doubleValue];
    double beatB = [beatMap[beatIndex + 1] doubleValue];
    // Main beat at zero, then fill in 3 subdivisions after the main beat.
    for (NSUInteger subdiv = 0; subdiv < 4; subdiv++) {
      [choppedBeatMap addObject:@( (beatA + (beatB - beatA) * (subdiv / 4.0)) )];
    }
  }
  // So that playback will play all of the last beat, add last beat and 3 subdivisions off the last
  // beat, extrapolating the deltas from the penultimate->last beat.
  double penultimateBeat = [beatMap[beatMap.count - 2] doubleValue];
  double lastBeat = [beatMap[beatMap.count - 1] doubleValue];
  double lastBeatSubdivDelta = (lastBeat - penultimateBeat) / 4.0;

  for (NSUInteger subdiv = 0; subdiv < 4; subdiv++) {
    [choppedBeatMap addObject:@( lastBeat + (subdiv * lastBeatSubdivDelta) )];
  }
  return choppedBeatMap;
}

- (BOOL)decompressAssetForURL:(NSURL *)url withBeatMap:(NSArray *)beatMap {
  NSArray *choppedBeatMap = [self chopBeatMap:beatMap]; //chop int 4 beats per beat

  // Turn into c-array of long sample values at playback rate.
  _beatMap = malloc(choppedBeatMap.count * sizeof(long));
  _beatMapCount = (int)choppedBeatMap.count;
  for (NSUInteger i = 0; i < _beatMapCount; i++) {
    _beatMap[i] = (long)([choppedBeatMap[i] doubleValue] * SAMPLE_RATE);
  }

  // Load asset.
  AVURLAsset *asset =
      [[AVURLAsset alloc] initWithURL:url
                              options:@{ AVURLAssetPreferPreciseDurationAndTimingKey:@(YES) }];
  Float64 durationSec = CMTimeGetSeconds([asset duration]); //sync load
  NSLog(@"Loaded asset of duration %f", durationSec);

  // Alloc audio data.
  if (_audioData) {
    free(_audioData);
  }
  _audioDataFrameCount = durationSec * SAMPLE_RATE;
  NSUInteger audioDataSampleCount = _audioDataFrameCount * NUM_CHANNELS;
  _audioData = (int16_t*)malloc(audioDataSampleCount * sizeof(int16_t));

  // Decompress audio into _audioData;
  __block NSUInteger audioDataIndex = 0;
  [BDPAudioUtils decompressAsset:asset
                    channelCount:NUM_CHANNELS
                      sampleRate:SAMPLE_RATE
                          handler:
   ^(int16_t *data, NSUInteger frameCount) {
    for (NSUInteger i = 0; i < frameCount * NUM_CHANNELS; i++, audioDataIndex++) {
      if (audioDataIndex >= audioDataSampleCount) {
        //NSLog(@"OFF THE END!");
        break;
      }
      _audioData[audioDataIndex] = data[i]; // HERE going off end of audio data for some reason...erratic? add check against
    }
   }];
  return YES;
}

#pragma mark -

- (void)beatDidChange {
  __weak AudioController *weakSelf = self;
  NSUInteger currBeat = _currBeat;
  dispatch_async(dispatch_get_main_queue(), ^{
    [weakSelf.delegate playbackHitBeatIndex:currBeat];
  });
}

- (void)playbackDidFinish {
  [self stop];
  __weak AudioController *weakSelf = self;
  dispatch_async(dispatch_get_main_queue(), ^{
    [weakSelf.delegate playbackDidFinish];
  });
}

@end
