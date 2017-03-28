#import "BDPAudioUtils.h"

#import <AVFoundation/AVFoundation.h>

@implementation BDPAudioUtils

+ (BOOL)decompressAsset:(AVURLAsset *)asset
           channelCount:(NSUInteger)channelCount
             sampleRate:(CGFloat)sampleRate
                handler:(void(^)(int16_t* data, NSUInteger frameCount))decompressHandler {
  NSError *error;
  AVAssetReader *myAssetReader = [[AVAssetReader alloc] initWithAsset:asset error:&error];
  if (error){
    NSLog(@"Error creating assetReader: %@",[error localizedDescription]);
  }

  NSDictionary *audioOptionsDictionary =
      @{ AVNumberOfChannelsKey : @(channelCount),
         AVFormatIDKey : @(kAudioFormatLinearPCM),
         AVLinearPCMBitDepthKey : @(16),
         AVLinearPCMIsFloatKey : @(NO),
         AVSampleRateKey : @(sampleRate) };

  NSArray<AVAssetTrack *> *audioTracks = [asset tracksWithMediaType:AVMediaTypeAudio];
  if (!audioTracks.count) {
    return NO; // No audio anywhere, so bad analaysis.
  }
  if (audioTracks.count > 1) {
    // If more than one, grab the first and hope...
    audioTracks = @[ audioTracks[0] ];
  }

  AVAssetReaderOutput *myOutput =
      [[AVAssetReaderAudioMixOutput alloc] initWithAudioTracks:audioTracks
                                                 audioSettings:audioOptionsDictionary];
  if (!myOutput){
    NSLog(@"Could not initialize the AVAssetReaderTrackOutput.");
    return NO;
  }

  if ([myAssetReader canAddOutput:myOutput]) {
    [myAssetReader addOutput:myOutput];
  } else {
    NSLog(@"Error: Cannot add output!!!");
    return NO;
  }

  if (![myAssetReader startReading]) {
    NSLog(@"Error: Asset reader cannot start reading. Error: %@",
          [myAssetReader.error localizedDescription]);
    return NO;
  }

  AudioBufferList buffList;
  int16_t* data;
  NSUInteger totalFramesRead = 0;
  CMBlockBufferRef blockBufferOut;

  while (YES) {
    CMSampleBufferRef myBuff = [myOutput copyNextSampleBuffer];
    if (!myBuff) {
      break;
    }
    NSUInteger frameCount = CMSampleBufferGetNumSamples(myBuff);
    totalFramesRead += frameCount;
    OSStatus err =
        CMSampleBufferGetAudioBufferListWithRetainedBlockBuffer(myBuff, NULL, &buffList,
            sizeof(buffList), NULL, NULL, 0, &blockBufferOut);

    if (err) {
      NSLog(@"CMSampleBufferGetAudioBufferListWithRetainedBlockBuffer error %d", (int)err);
      if(blockBufferOut) {
        CFRelease(blockBufferOut);
      }
      return NO;
    }

    data = (int16_t*)buffList.mBuffers[0].mData;
    if(myBuff) {
      CFRelease(myBuff);
    }
    if(frameCount == 0) {
      NSLog(@"read %lu samples", (unsigned long)totalFramesRead);
      break;
    }

    if (decompressHandler) {
      decompressHandler(data, frameCount);
    }

    if(blockBufferOut) {
      CFRelease(blockBufferOut);
    }
  }
  return YES;
}

@end
