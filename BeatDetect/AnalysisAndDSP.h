#include <Accelerate/Accelerate.h>

// Initialize.
void AccelInit(int analysisFrameSize, int sampleRate, int hopSize, int melSpecSize);
// Deallocate.
void AccelDealloc();
// Return melodic spectrum of the input data.
// input is expected to be of length analysisFrameSize. output is of size melSpecSize.
float* AccelRender(short* ioData);

/// Analyze an onset function for beat information. Returns int array of frame offsets of the
/// detected beats. Populates outputLength and ratio DEI
//long* checkOnsets(double* onsetFunction, long onsetLength, long *outputLength, int* subdivision);

/// Take the output of checkOnsets (as inputData) as well as the onset function values at those
/// beats (as onsetVals), and do additional processing to "fold" the BPM (and all the beat values)
/// between a desired range of values. Currently folds between 100 and 200 BPM.
//long* massageBeats(long* inputData, long inputDataLength, int subdivision, float* onsetVals,
  //                 long* outputDataLength);

long* getBeatMap(double* onsetFunction, long onsetFunctionLength, long* outputLength);
