#import <Foundation/Foundation.h>

@interface AnalysisController : NSObject

// Return a beat map of NSNumber floats, representing position, in seconds, of each beat.
- (NSArray<NSNumber *> *)analyzeAssetForURL:(NSURL *)url;

@end
