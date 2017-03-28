#define BPM_COMPARE_EPSILON 1.0
#define BEATCOUNT_COMPARE_EPSILON 4

#import <XCTest/XCTest.h>

#import "AnalysisController.h"
#import "AudioController.h"

@interface BeatDetectTests : XCTestCase
@end


@implementation BeatDetectTests

- (void)testGenerateBaseline {
  // Generate info for testing.
  NSArray *filePaths =
      [[NSBundle bundleForClass:[self class]] pathsForResourcesOfType:@".mp3" inDirectory:nil];
  for(NSUInteger i = 0; i < [filePaths count]; i++) {
    NSString *path = [filePaths objectAtIndex:i];

    NSArray *beatMap = [self beatMapForPath:path];
    CGFloat avgBPM = [self avgBPMForBeatMap:beatMap];
    NSString *key = [path lastPathComponent];
    printf("\nexpectedBeatCountMap[@\"%s\"] = @(%d);",
        [key cStringUsingEncoding:NSUTF8StringEncoding], (int)[beatMap count]);
    printf("\nexpectedBpmMap[@\"%s\"] = @(%.3f);",
        [key cStringUsingEncoding:NSUTF8StringEncoding], (float)avgBPM);
  }
}

- (void)testGenerateTestDataMap {
  NSMutableDictionary *expectedBeatCountMap = [NSMutableDictionary dictionary];
  NSMutableDictionary *expectedBpmMap = [NSMutableDictionary dictionary];

  // Add expected output for test audio files in bundle.
  //expectedBeatCountMap[@"filename.mp3"] = @(750);
  //expectedBpmMap[@"filename.mp3"] = @(156.050);

  NSArray *filePaths =
      [[NSBundle bundleForClass:[self class]] pathsForResourcesOfType:@".mp3" inDirectory:nil];

  NSUInteger beatCountPassCount = 0;
  NSUInteger bpmPassCount = 0;
  NSUInteger fileCount = [filePaths count];
  for(NSUInteger i = 0; i < fileCount; i++) {
    NSString *path = [filePaths objectAtIndex:i];
    NSString *key = [path lastPathComponent];
    NSLog(@"testing %@", key);
    XCTAssertNotNil(expectedBeatCountMap[key]);
    XCTAssertNotNil(expectedBpmMap[key]);

    NSArray *beatMap = [self beatMapForPath:path];
    CGFloat avgBPM = [self avgBPMForBeatMap:beatMap];

    NSLog(@"expected beat count %lu bpm %.2f: actual %lu, %.2f",
          [expectedBeatCountMap[key] integerValue], [expectedBpmMap[key] floatValue],
          (unsigned long)beatMap.count, avgBPM);

    BOOL beatCountIsEqual =
        labs((NSInteger)[beatMap count] - [expectedBeatCountMap[key] integerValue]) <
        BEATCOUNT_COMPARE_EPSILON;
    XCTAssertTrue(beatCountIsEqual);
    if (beatCountIsEqual) {
      beatCountPassCount++;
    }
    BOOL bpmIsEqual = fabs(avgBPM - [expectedBpmMap[key] floatValue]) < BPM_COMPARE_EPSILON;
    XCTAssertTrue(bpmIsEqual);
    if (bpmIsEqual) {
      bpmPassCount++;
    }
    /*if (avgBPM > 0) {
      [self play:path beatMap:beatMap];
    }*/
  }
  NSLog(@"====beat count pass %lu / %lu, bpm pass %lu / %lu", (unsigned long)beatCountPassCount,
        (unsigned long)fileCount, (unsigned long)bpmPassCount, (unsigned long)fileCount);
}

- (void)play:(NSString *)path beatMap:(NSArray *)beatMap beatOffset:(NSUInteger)beatOffset {
  XCTAssertNotNil(path);
  AudioController *ac = [[AudioController alloc] init];
  NSURL *url = [NSURL fileURLWithPath:path];
  ac.playbackClicks = YES;
  [ac decompressAssetForURL:url withBeatMap:beatMap];
  [ac setBPM:0];
  [ac jumpToBeatIndex:beatOffset];
  [ac play];
  [NSThread sleepForTimeInterval:20];
  [ac stop];
}

- (NSArray *)beatMapForPath:(NSString *)path {
  NSURL *url = [NSURL fileURLWithPath:path];
  AnalysisController *ac = [[AnalysisController alloc] init];
  return [ac analyzeAssetForURL:url];
}

- (CGFloat)avgBPMForBeatMap:(NSArray *)beatMap {
  if ([beatMap count] < 2) {
    return 0;
  }

  // diff between first beat and last beat. First beat can have artifacts, so skip that.
  double deltaSum = [beatMap[beatMap.count - 1] doubleValue] - [beatMap[0] doubleValue];
  // Get average beat length.
  // Minus one because we are not taking the full duration of the last beat into account
  // (since we only know its start, not its end).
  deltaSum /= ([beatMap count] - 1);

  return 60.0/deltaSum;
}

@end
