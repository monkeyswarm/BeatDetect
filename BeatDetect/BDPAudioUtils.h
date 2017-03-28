#import <Foundation/Foundation.h>

@class AVURLAsset;

@interface BDPAudioUtils : NSObject

/**
 * Decompress AVURLAsset into raw audio frames. Decompresses to 16 bit integer linear PCM, with
 * provided channel count and sample rate.
 * Handler is called synchronously while there are chunks of decompressed audio. |data| is of size
 * frameCount * channelCount.
 */
+ (BOOL)decompressAsset:(AVURLAsset *)asset
           channelCount:(NSUInteger)channelCount
             sampleRate:(CGFloat)sampleRate
                handler:(void(^)(int16_t* data, NSUInteger frameCount))decompressHandler;

@end
