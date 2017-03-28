#import "AnalysisController.h"

#import <AVFoundation/AVFoundation.h>

#import "AnalysisAndDSP.h"
#import "BDPAudioUtils.h"

#define ANALYSIS_SR 44100
#define ANALYSIS_NPOINTS 1024
#define HOP_SIZE 512
#define melSpecDim 40

@implementation AnalysisController

- (instancetype)init {
  self = [super init];
  if (self) {
    AccelInit(ANALYSIS_NPOINTS, ANALYSIS_SR, HOP_SIZE, melSpecDim);
  }
  return self;
}

- (void)dealloc{
  AccelDealloc();
}

- (NSArray<NSNumber *> *)analyzeAssetForURL:(NSURL *)url {
  AVURLAsset *asset =
      [[AVURLAsset alloc] initWithURL:url
                              options:@{ AVURLAssetPreferPreciseDurationAndTimingKey:@(YES) }];
  Float64 durationSec = CMTimeGetSeconds([asset duration]); // Loads asset sync'ly.

  // Analysis data
  __block NSUInteger frameIndex = 0;
  __block NSUInteger shuttleIndex = 0;
  float *prevSpec = (float*)malloc(melSpecDim * sizeof(float));
  memset(prevSpec, 0, melSpecDim*sizeof(float)); //neccesary?

  NSUInteger estimatedFrames = (durationSec + 1) * ANALYSIS_SR / HOP_SIZE;
  double* onsetFunction = (double *)malloc(estimatedFrames * sizeof(double));
  memset(onsetFunction, 0, estimatedFrames * sizeof(double));

  int16_t* shuttle = (int16_t *)malloc(ANALYSIS_NPOINTS * sizeof(int16_t));
  memset(shuttle, 0, ANALYSIS_NPOINTS * sizeof(int16_t));

  [BDPAudioUtils decompressAsset:asset
                    channelCount:1
                      sampleRate:ANALYSIS_SR
                         handler:
    ^(int16_t *data, NSUInteger frameCount) {
      for (NSUInteger i = 0; i < frameCount; i++, shuttleIndex++) {
        shuttle[shuttleIndex] = data[i];
        if (shuttleIndex == ANALYSIS_NPOINTS - 1) { // We just put in last sample into shuttle.
          // Get the melodic spectrum of the contents of the shuttle.
          // Note that currSpec is owned (and reused) by AccelRender, not allocated each call.
          float* currSpec = AccelRender(shuttle);
          for (int j = 0; j < melSpecDim; j++) {
            if (frameIndex > 0) {
              onsetFunction[frameIndex] += fabs(currSpec[j] - prevSpec[j]);
            }
            prevSpec[j] = currSpec[j];
          }
          onsetFunction[frameIndex] /= melSpecDim;
          for (int j = HOP_SIZE; j < ANALYSIS_NPOINTS; j++) {
            shuttle[j - HOP_SIZE] = shuttle[j];
          }
          // Reset the position of shuttleIndex to accept upcomgin data.
          shuttleIndex = (ANALYSIS_NPOINTS-HOP_SIZE) - 1; //-1! cause we increment it next time!
          frameIndex++;
        }
      }
    }];

  long beatDataLength = 0;
  long* beatData = getBeatMap(onsetFunction, frameIndex, &beatDataLength);

  if(beatData == nil || beatDataLength == 0) {
    printf("\nCould not get beats from getBeatMap.");
    return nil;
  }

  NSMutableArray* beatMap = [[NSMutableArray alloc] initWithCapacity:beatDataLength];

  // Convert massaged data (in fft frames) to time in seconds.
  double multiplier = (double)HOP_SIZE / ANALYSIS_SR;
  for(NSUInteger i = 0; i < beatDataLength; i++) {
    [beatMap addObject:[NSNumber numberWithDouble:(beatData[i] * multiplier)]];
  }

  free(onsetFunction);
  free(shuttle);
  free(prevSpec);
  free(beatData);
  return beatMap;
}

@end
