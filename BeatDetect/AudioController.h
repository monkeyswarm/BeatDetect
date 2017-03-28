#import <Foundation/Foundation.h>


///
@protocol AudioControllerDelegate <NSObject>

- (void)playbackHitBeatIndex:(NSUInteger)beatIndex;

- (void)playbackDidFinish;

@end

///
@interface AudioController : NSObject

@property(nonatomic, weak) id<AudioControllerDelegate> delegate;
@property(nonatomic, getter=hasPlaybackClicks) BOOL playbackClicks;

- (BOOL)decompressAssetForURL:(NSURL *)url withBeatMap:(NSArray *)beatMap;
- (void)play;
- (void)stop;
- (void)setBPM:(CGFloat)BPM; // Set to zero to play back at normal speed.
- (void)jumpToBeatIndex:(NSUInteger)beatIndex;

@end
