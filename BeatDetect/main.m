#import <Foundation/Foundation.h>

#import "AnalysisController.h"
#import "AudioController.h"

@interface AppDelegate : NSObject <AudioControllerDelegate>
@end

@implementation AppDelegate {
  NSUInteger _beatCount;
}

- (instancetype)init {
  self = [super init];

  return self;
}

- (void)showUsage {
  printf("\nUsage:\n    -file <relative file path>\n    -shouldPlay <YES or NO> (defaults to YES)\n    -BPM <float> (playback BPM value, defaults to 0, which plays back at original speed)\n    -shouldPrintToScreen <YES or NO> (whether to print beat map to stdout)\n    -outputFile <relative file path> (output file to write beat map)\n    -shouldClick <YES or NO> (whether to add a click to playback to hear beat analysis)\n    -startBeat <integer> (beat index on which to start playback)\n");
}

- (int)run {
  // Get command line arguments from standard defaults.
  NSUserDefaults *standardDefaults = [NSUserDefaults standardUserDefaults];

  // check for nothing to do.
  if (![standardDefaults boolForKey:@"shouldPrintToScreen"] &&
      ![standardDefaults boolForKey:@"outputFile"] &&
      ([standardDefaults objectForKey:@"shouldPlay"] && ![standardDefaults boolForKey:@"play"])) {
    printf("\nYou aren't doing anything with the analysis (printing to screen, writing to file, or playing audio). Exiting.");
    [self showUsage];
    return 1;
  }

  //file
  NSString *relativeFilePath = [standardDefaults stringForKey:@"file"];
  relativeFilePath =
      [relativeFilePath stringByAddingPercentEscapesUsingEncoding:NSUTF8StringEncoding];
  if (!relativeFilePath) {
    [self showUsage];
    return 1;
  }

  NSString *currentPath = [[NSFileManager defaultManager] currentDirectoryPath];
  NSURL *currentPathFileURL = [NSURL fileURLWithPath:currentPath];

  NSURL *url = [NSURL URLWithString:relativeFilePath relativeToURL:currentPathFileURL];
  if (!url) {
    NSLog(@"Could not generate file url. \n %@ \n %@", relativeFilePath, currentPathFileURL);
    return 1;
  }
  if (![[NSFileManager defaultManager] fileExistsAtPath:[url path]]) {
    NSLog(@"Could not find file at path:%@", [url path]);
    return 1;
  }

  // Analyze.
  AnalysisController *ac = [[AnalysisController alloc] init];
  NSArray *beatMap = [ac analyzeAssetForURL:url];
  _beatCount = beatMap.count;
  CGFloat avgNativeBpm = [self avgBPMForBeatMap:beatMap];
  printf("\nComputed beat map of %lu beats. Average native BPM = %.2f",
         (unsigned long)beatMap.count, avgNativeBpm);

  // Handle data output.
  BOOL shouldPrintOutput = [standardDefaults boolForKey:@"shouldPrintToScreen"];
  NSString *outputRelativeFilePath = [standardDefaults stringForKey:@"outputFile"];
  outputRelativeFilePath =
      [outputRelativeFilePath stringByAddingPercentEscapesUsingEncoding:NSUTF8StringEncoding];
  NSURL *outputURL = [NSURL URLWithString:outputRelativeFilePath relativeToURL:currentPathFileURL];
  NSString *outputFilePath = [outputURL path];

  if (shouldPrintOutput) {
    printf("\nbeat map:");
    for (NSNumber *beatOffset in beatMap) {
      printf("\n%f", [beatOffset floatValue]);
    }
    printf("\n");
  }

  if (outputRelativeFilePath) {
    NSMutableString *outputString = [NSMutableString string];
    for (NSNumber *beatOffset in beatMap) {
      [outputString appendString:[NSString stringWithFormat:@"%f\n", [beatOffset floatValue]]];
    }
    NSError *error = nil;
    [outputString writeToFile:outputFilePath
                   atomically:YES
                     encoding:NSUTF8StringEncoding
                        error:&error];
    if (error) {
      printf("\nCould not write to output file path %s",
             [outputFilePath cStringUsingEncoding:NSUTF8StringEncoding]);
      return 1;
    }
  }

  if ([beatMap count] < 2) {
    NSLog(@"Could not detect beat content for this audio.");
    return 1;
  }

  // playback
  BOOL shouldPlay = YES; // Default to YES
  if ([standardDefaults objectForKey:@"shouldPlay"]) {
    shouldPlay = [standardDefaults boolForKey:@"shouldPlay"];
  }

  if (shouldPlay) {
    BOOL shouldClick = [standardDefaults boolForKey:@"shouldClick"];
    CGFloat userBPM = [standardDefaults floatForKey:@"BPM"];
    NSInteger userStartBeat = [standardDefaults integerForKey:@"startBeat"];
    AudioController *audioController = [[AudioController alloc] init];
    audioController.delegate = self;

    BOOL decompressSuccess = [audioController decompressAssetForURL:url withBeatMap:beatMap];
    if (!decompressSuccess) {
      return 1;
    }
    [audioController setBPM:userBPM]; // 0 is native speed.
    if (userStartBeat > 0) {
      [audioController jumpToBeatIndex:userStartBeat];
    }
    audioController.playbackClicks = shouldClick;

    if (userBPM > 0) {
      printf("\nplaying back at BPM %.2f", userBPM);
    } else {
      printf("\nplaying back at native speed");
    }

    printf("\n"); // New line so beat printing is on new line.

    [audioController play];
    // Keep program alive while audio plays. Main thread will wait until |playbackDidFinish|.
    [[NSRunLoop currentRunLoop] run];
  }

  // When shouldPlay is YES, this is never reached. Finishing audio exits in |playbackDidFinish|.
  return 0;
}

#pragma mark - AudioControllerDelegate

- (void)playbackHitBeatIndex:(NSUInteger)beatIndex {
  printf("beat %lu / %lu \r", (unsigned long)beatIndex, (unsigned long)_beatCount);
  fflush(__stdoutp);
}

- (void)playbackDidFinish {
  printf("\n");
  CFRunLoopStop(CFRunLoopGetCurrent());
  // Stopping the run loop doesn't continue the [run] method to completion. Force exit.
  exit(0);
}

#pragma mark - Private

- (CGFloat)avgBPMForBeatMap:(NSArray *)beatMap {
  if ([beatMap count] < 2) {
    return 0;
  }

  // Difference between first and last beat.
  double deltaSum = [beatMap[beatMap.count - 1] doubleValue] - [beatMap[0] doubleValue];
  // Get average beat length.
  // Minus one because we are not taking the full duration of the last beat into account
  // (since we only know its start, not its end).
  deltaSum /= (beatMap.count - 1);

  return 60.0 / deltaSum;
}

@end


int main(int argc, const char * argv[]) {
  int result = 0;
  @autoreleasepool {
    AppDelegate *appDelegate = [[AppDelegate alloc] init];
    result = [appDelegate run];
  }
  return result;
}
