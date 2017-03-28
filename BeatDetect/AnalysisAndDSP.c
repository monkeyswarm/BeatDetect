#include "AnalysisAndDSP.h"

//
//#include <float.h>
//#include <math.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <limits.h>

#define minBPM 100
#define maxBPM 200 // This should be 2x the minBPM, otherwise unusual folding may occur.
// TODO: future work will allow arbitrary BPM ranges.

static const bool verbose = false;

static const float rmsTarget = .1; // Desired RMS that the volume adjustment will target.
static const float rmsAlpha = 0.001; // Smoothing factor for RMS volume adjustment.
static const float melSpecSilenceThresh = -40.0; // vol threshold (in dB) for a melSpec band to be silent.

// Use to define the minimum frequency bin to start scanning for initial xcorr peak. Converted to a
// bin index based on the sampling rate and hop size. For 44100 & 512, this should become bin 10.
static const int onsetSmallestPeakThresholdBPM = 516;

// Use to define the minimum frequency bin, below which we MUST find a valid duple/triple/quadruple,
// else it is not a valid analysis. For 44100 & 512, this should become bin 20.
static const int onsetSmallestMultThresholdBPM = 258;

// Use to define the frequency bin up to which we look for tuples. For 44100 & 512, this should
// become bin 30. TODO: make this tuple-specific (e.g. bpm thresh for duples vs triplets).
//static const int onsetLookForTuplesMaxBPM = 172;

// The reported periodicity peak must have a value greater than this to be reported as valid.
static const int minValidPeakVal = 1200;
static const int minValidAdjustedPeakVal = 600;

// Number of bins to scan in the autocorrelation
static const int autoCorrelationMaxLag = 128;

// FFT and output info.
int _analysisFrameSize;
int _sampleRate;
int _hopSize;
int _melSpecSize;
FFTSetup _fftSetup;

// Struct holding real/imag components of a single FFT analysis.
DSPSplitComplex _complexComponents;

// Holds moving rms value. Starts at rmsTarget.
float _rms;
// Holds a frame of time-domain samples for pre-FFT processing.
float* _floatShuttle;
// Holds the N-band melodic spectrum that will be returned from each frame.
float* _melSpec;
// A double version of the above for use pre-decibel computation.
double* _melSpecD;
// Array
double *_melOfLin;//will be array of N holding mappings to MelSpec
double *_melCenter;//hold 40
double *_melWidth;//hold 40

// pointer to an array of size _melSpecDim, of arrays of weights per bin.
double** _ptrArrayToWeightArrays;
int* _weightArrayCount;
int* _weightArrayStartBin;

float* _window;

// Library of DSP functions. Adapted from DSP.java in meapsoft.
// https://www.ee.columbia.edu/~ronw/code/MEAPsoft/src/com/meapsoft/DSP.java

// Convolves sequences a and b.  The resulting convolution has length a.length+b.length-1.
double* conv(double* a, long alength, double* b, long blength) {
	long ylength = alength + blength - 1;
	double* y = (double*)malloc(ylength * sizeof(double));

	// make sure that a is the shorter sequence
	if (alength > blength) {
    // Swap values.
    double* tmp = a;
		a = b;
		b = tmp;
    long templen;
		templen = alength;
		alength = blength;
		blength = templen;
	}

	for (int lag = 0; lag < ylength; lag++) {
		y[lag] = 0;
		// where do the two signals overlap?
		long start = 0;
		// we can't go past the left end of (time reversed) a
    if (lag > alength - 1) {
			start = lag - alength + 1;
    }
		long end = lag;
		// we can't go past the right end of b
    if (end > blength - 1) {
			end = blength - 1;
    }
		for (long n = start; n <= end; n++) {
			y[lag] += b[n] * a[lag - n];
		}
	}
	
	return y;
}

// cross correlation of a and b
double* xcorr(double *a, long alength, double *b, long blength, int maxlag) {
	int ylength = 2 * maxlag + 1;
	double *y = (double*)malloc(ylength * sizeof(double));
  memset(y, 0, ylength * sizeof(double));

	for (long lag = blength - 1, idx = maxlag - blength + 1; lag > -alength; lag--, idx++) {
    if (idx < 0) {
			continue;
    }
		
    if (idx >= ylength) {
			break;
    }
		
		// where do the two signals overlap?
		long start = 0;
		// we can't start past the left end of b
		if (lag < 0) {
			start = -lag;
		}

		long end = alength - 1;
		// we can't go past the right end of b
		if (end > blength - lag - 1) {
			end = blength - lag - 1;
		}
		
		for (long i = start; i <= end; i++) {
			y[idx] += a[i] * b[lag + i];
		}
	}
	
	return y;
}

// Filters x by the IIR filter defined by a and b.
double* bdp_filter(double* b, long blength, double* a, long alength, double* x, long xlength) {
	double* y = (double*)malloc(xlength * sizeof(double));
  memset(y, 0, xlength * sizeof(double));

	// factor out a[0]
	if (a[0] != 1) {
		for (int ia = 1; ia < alength; ia++)
			a[ia] = a[ia] / a[0];
		
		for (int ib = 0; ib < blength; ib++)
			b[ib] = b[ib] / a[0];
	}

	for (int t = 0; t < xlength; t++) {
		y[t] = 0;
		
		// input terms
		long len = blength - 1 < t ? blength - 1 : t;
		for (int ib = 0; ib <= len; ib++) {
			y[t] += b[ib] * x[t - ib];
		}
		// output terms
		len = alength - 1 < t ? alength - 1 : t;
		for (int ia = 1; ia <= len; ia++) {
			y[t] -= a[ia] * y[t - ia];
		}
	}
	return y;
}

// Sum of the contents of a.
double sum(double* a, int alength) {
	double y = 0;
  for (int x = 0; x < alength; x++) {
		y += a[x];
  }
	return y;
}

// Multiplies each element of a[] by b.
double* bdp_times_array_scalar(double *a, int alength, double b) {
	double *y = (double *)malloc(alength*sizeof(double));
  for (int x = 0; x < alength; x++) {
		y[x] = a[x] * b;
  }
	return y;
}

// Element-wise mult of a and b. a and b must be of the same length.
double* bdp_times_array_array(double* a, int alength, double* b) {
	double* y = (double*)malloc(alength * sizeof(double));
  for (int x = 0; x < alength; x++) {
		y[x] = a[x]*b[x];
  }
	return y;
}

// Returns the max element of a.
double bdp_max(double* a, long alength) {
  double y = -DBL_MAX;
  for (long x = 0; x < alength; x++) {
    if (a[x] > y) {
      y = a[x];
    }
  }
  return y;
}

// Returns the index of the max element of a.
int argmax(double* a, int alength) {
	double y = -DBL_MAX;
	int idx = -1;
	int x;
	
	for (x = 0; x < alength; x++) {
		if (a[x] > y) {
			y = a[x];
			idx = x;
		}
	}
	
	return idx;
}

// Returns (new copy) of the slice of array a between start and end, inclusive.
// start/end must be valid, there is no length checking.
double* slice(double* a, long start, long end) {
	double* y = (double*)malloc((end - start + 1) * sizeof(double));

  for (long x = start, iy = 0; x <= end; x++, iy++) {
		y[iy] = a[x];
  }
	
	return y;
}

// Returns an array as follows: {start, start+increment, ... , end-increment, end}
double* range(int start, int end, int increment) {
	int ylength = (1 + (end - start) / increment);
	double* y = (double*)malloc(ylength * sizeof(double));
  for (int x = 0, num = start; x < ylength; x++, num += increment) {
		y[x] = num;
  }
	return y;
}

// Return array of values from a, from the indeces of idx the value of idx was 1.
double* subsref(double* a, long alength, int* idx, long idxlength, int *subsreflength) {
	int ones_count = 0;
  // Determine length of results by counting 1s in idx.
  for (long x = 0; x < idxlength; x++) {
    if (idx[x] > 0) {
      ones_count++;
    }
  }

  double* y = (double*)malloc(ones_count * sizeof(double));

  int y_idx = 0;
  for (int x = 0; x < idxlength; x++) {
    if (idx[x] > 0) {
      y[y_idx++] = a[x]; // get value of index in a. assumes alength==idxlength.
    }
  }

  *subsreflength = ones_count;
	return y;
}

// Sorting function for qsort
int compare_doubles (const void *a, const void *b) {
	double temp = *(const double*)a - *(const double*)b; // cast to double.
  if (temp > 0) {
		return 1;
  } else if (temp < 0) {
		return -1;
  } else {
		return 0;
  }
}

// Returns the lower median value contained in a.
double median(double* a, int alength) {
  // Make a copy to sort.
	double* tmp = slice(a, 0, alength - 1);
	qsort(tmp, alength, sizeof (double), compare_doubles);
	double toReturn = tmp[(int)alength / 2];
	free(tmp);
	return toReturn;
}

// Raise each element of a to the b'th power.
double* power(double* a, int alength, double b) {
	double* y = (double *)malloc(alength * sizeof(double));
  for (int x = 0; x < alength; x++) {
		y[x] = pow(a[x], b);
  }

	return y;
}

// Take the natural log of each element of a.
double* bdp_log(double* a, int alength) {
	double* y = (double *)malloc(alength * sizeof(double));
  for (int x = 0; x < alength; x++) {
		y[x] = log(a[x]);
  }
	
	return y;
}

// Take the natural exp of each element of a.
double* bdp_exp(double* a, int alength) {
	double* y = (double *)malloc(alength * sizeof(double));
  for (int x = 0; x < alength; x++) {
		y[x] = exp(a[x]);
  }
	
	return y;
}

// Find local maxima in function a. Returns an array with 1's where a has a local maximum.
// Result is of the same length as input.
// TODO consider using byte/char/boolean type.
int* localmax(double* a, long alength) {
	int* v =(int*)malloc(alength * sizeof(int));
  memset(v, 0, alength * sizeof(int));

	for (long fr = 1; fr < alength - 1; fr++) {
		if (a[fr] > (a[fr - 1]) && a[fr] > (a[fr + 1])) {
			v[fr] = 1;
		}
	}
	
	return v;
}

// Convert a BPM into a bin index in the peaks array.
int binIndexForBPM(int BPM) {
  // E.g for sample rate 44100 and hop size 512:
  // bin = (60 / minBPM) * (44100 / 512), which turns into (44100 * 60) / (minBPM * 512).
  return (int)((_sampleRate * 60) / ((float)BPM * _hopSize));
}

// Take the onset function and return a list of beat positions.
// Adapted from Dan Ellis' beat detection.
// https://www.ee.columbia.edu/~ronw/code/MEAPsoft/src/com/meapsoft/DpweBeatOnsetDetector.java
// "Beat Tracking by Dynamic Programming"
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.94.8773&rep=rep1&type=pdf
// This adaptation adds additional computation to compute duple/triplet/quadruplets for the base
// tempo, rather than use the smallest significant peak (i.e. a faster BPM, which may be a
// subdivision of the beat). These results are additionally

long* getBeatMap(double* onsetFunction, long onsetFunctionLength, long* outputLength) {
  int subdivision = 0;
  // Remove DC
  double b[] = { 1, -1 };
  double a[] = { 1, -0.99 };
	double* fmm = bdp_filter(b, 2, a, 2, onsetFunction, onsetFunctionLength);

  // Autocorrelation.
	double* xfmm2 = xcorr(fmm, onsetFunctionLength, fmm, onsetFunctionLength, autoCorrelationMaxLag);

	// Find local maxima in the global autocorrelation (from the second half of the results).
	int xfmmlength = autoCorrelationMaxLag + 1; // Resulting size of slice().

  // Take second half of xcorr.
	double* xfmm = slice(xfmm2, autoCorrelationMaxLag, 2 * autoCorrelationMaxLag);
  free(xfmm2);
	int* xpks = localmax(xfmm, xfmmlength); // Fill xpks with value 1 at local maxima.

	// xpks should not include 'edge peak' at first index, but ensure it.
	xpks[0] = 0;

	// Largest local max after first neg pts.
  int subsreflength = 0;
	double* temp_subref = subsref(xfmm, xfmmlength, xpks, xfmmlength, &subsreflength);
  double maxPeakValue = bdp_max(temp_subref, subsreflength);
  free(temp_subref);

  if (maxPeakValue < minValidPeakVal) {
    if (verbose) {
      printf("\nmax peak not strong enough!");
    }
    *outputLength = 0;
    return NULL;// thisbeat;
  }
	
	// Then period is shortest period with a peak that approaches the max. It shouldn't be too small.
  // Find first peak of >60% of max peak.
  int pd = 0;
  int startLookIndex = binIndexForBPM(onsetSmallestPeakThresholdBPM);
  for (int x = startLookIndex; x < xfmmlength; x++) {
		if (xpks[x] == 1) {
			if (xfmm[x] > 0.6 * maxPeakValue) {
				pd = x;
				break;
			}
		}
	}

  // look at peaks
  if (verbose) {
    printf("\npks:");
    for (int x = 0; x < xfmmlength; x++) {
      if (xpks[x] == 1) {
        printf(" %d[%.2f]", x, xfmm[x]);
      }
    }
    printf("\nPERIOD:%d", pd);
  }

  if (pd == 0) {
    if (verbose) {
      printf("failed to find good peak after %d", startLookIndex);
    }
    *outputLength = 0;
    return NULL;
  }

	// In addition to the peak we've found, if it is fairly small period (fast BPM), look at the peaks
  // of the autocorrelation and see if there is a nice peak at a multiple of this one.
  // Look at double, triple, and quadruple peaks. If used, store multiplier as subdivision.
  // OR if the period is too large (slow BPM), compute the subdivision (but keep the slow period).
  subdivision = 1;
  //int lookForTuplesUnderBin = binIndexForBPM(onsetLookForTuplesMaxBPM);
  int maxBPMBin = binIndexForBPM(maxBPM);
  int minBPMBin = binIndexForBPM(minBPM) + 1;

  if (pd <= maxBPMBin) { // period is too fast, multiply
    // For duple, triple, and quad, see if there is a peak at 2,3, and 4 times the base peak.
    double dupleVal = 0;
    int dupleIndex = -1;
    for (int i = (pd - 1) * 2; i <= (pd + 1) * 2; i++) {
      if (xpks[i] == 1) {
        dupleVal = xfmm[i];
        dupleIndex = i;
      }
    }

    double tripleVal = 0;
    int tripleIndex = -1;
    for (int i = (pd - 1) * 3; i <= (pd + 1) * 3; i++) {
      if (xpks[i] == 1) {
        tripleVal = xfmm[i];
        tripleIndex = i;
      }
    }

    double quadVal = 0;
    int quadIndex = -1;
    for (int i = (pd - 1) * 4; i <= (pd + 1) * 4; i++) {
      if (xpks[i] == 1) {
        quadVal = xfmm[i];
        quadIndex = i;
      }
    }

    // Look for, in this order: duple, triple, and quadruple.
    // If duple peak is also >50% of max peak, use it.
    if ((dupleIndex > -1) && (xfmm[dupleIndex] > .5 * maxPeakValue)) {
      pd = dupleIndex;
      subdivision = 2;
    }

    // If there is a duple peak and triple peak, but the triple peak is much stronger, use it.
    // Or, if there is not a duple peak but the triple peak is >60% of max peak, use it.
    if (((tripleIndex > -1 && dupleVal == -1) || ((float)tripleVal / dupleVal > 1.5)) &&
        xfmm[tripleIndex] > .6 * maxPeakValue) {
      pd = tripleIndex;
      subdivision = 3; // So playback knows to triple speed instead of double.
    }

    // If there is a quadruple peak, and it is greater than triple val, and it is >50% of max peak,
    // use it.
    if ((quadIndex > -1) && (xfmm[quadIndex] > .5 * maxPeakValue) && (quadVal > tripleVal)) {
      pd = quadIndex;
      subdivision = 4;
    }
  } else if (pd >= minBPMBin) {
    // it's too slow, don't change the period (keep it slow) but record the subdivision.
    // For duple, triple, and quad, see if there is a peak at 2,3, and 4 fraction of the base peak.
    double dupleVal = 0;
    int dupleIndex = -1;
    for (int i = (int)floor((pd - 1) / 2.0); i <= (int)ceil((pd + 1) / 2.0); i++) {
      if (xpks[i] == 1) {
        dupleVal = xfmm[i];
        dupleIndex = i;
      }
    }

    double tripleVal = 0;
    int tripleIndex = -1;
    for (int i = (int)floor((pd - 1) / 3.0); (int)ceil(i <= (pd + 1) / 3.0); i++) {
      if (xpks[i] == 1) {
        tripleVal = xfmm[i];
        tripleIndex = i;
      }
    }

    double quadVal = 0;
    int quadIndex = -1;
    for (int i = (int)floor((pd - 1) / 4.0); (int)ceil(i <= (pd + 1) / 4.0); i++) {
      if (xpks[i] == 1) {
        quadVal = xfmm[i];
        quadIndex = i;
      }
    }

    // Look for, in this order: duple, triple, and quadruple.
    // If duple peak is also >30% of max peak, use it.
    // (Note, these percents are more permissive than when folding in the other direction above).
    if ((dupleIndex > -1) && (xfmm[dupleIndex] > .3 * maxPeakValue)) {
      subdivision = 2; // So massageBeats will look for 2x speed.
    }
    // If there is a duple peak and triple peak, but the triple peak is much stronger, use it.
    // Or, if there is not a duple peak but the triple peak is >30% of max peak, use it.
    if (((tripleIndex > -1 && dupleVal == -1) || ((float)tripleVal / dupleVal > 1.5)) &&
        xfmm[tripleIndex] > .3 * maxPeakValue) {
      subdivision = 3; // So massageBeats will look for 3x speed.
    }
    // If there is a quadruple peak, and it is greater than triple val, and it is >30% of max peak,
    // use it.
    if ((quadIndex > -1) && (xfmm[quadIndex] > .3 * maxPeakValue) && (quadVal > tripleVal)) {
      subdivision = 4; // So massageBeats will look for 4x (or 2x) speed.
    }
  }

  // Check that our new tempo peak is above our strength threshold.
  if (xfmm[pd] < minValidAdjustedPeakVal) {
    if (verbose) {
      printf("\npd peak not strong enough!");
    }
    *outputLength = 0;
    return NULL;
  }

  // Check that our new tempo peak is not too small (i.e. very high BPM and couldn't find a tuple).
  if (pd < binIndexForBPM(onsetSmallestMultThresholdBPM)) {
    if (verbose) {
      printf("fail :small tick with no valid tuple peak!");
    }
    *outputLength = 0;
    return NULL;
  }

  if (verbose) {
    printf("\n===Subdivision: %d beat peak %d (%.2f) max peak val (%.2f)\n",
           subdivision, pd, xfmm[pd], maxPeakValue);
  }
	
	// Smooth beat events
	int templtlength = (1 + (pd - -pd)); //after "range()", all others should keep length the same
	double* temp_range = range(-pd, pd, 1); // create -pd to pd ramp.
	double* temp_times = bdp_times_array_scalar(temp_range, templtlength, 32.0 / pd);
	double* temp_power = power(temp_times , templtlength, 2.0);
	double* temp_times2 = bdp_times_array_scalar(temp_power, templtlength, -0.5);
	double* templt = bdp_exp(temp_times2, templtlength);
	free(temp_range);
	free(temp_times);
	free(temp_power);
	free(temp_times2);
	
	long localscorelength = templtlength + onsetFunctionLength + 1;
	double* localscore_pre = conv(templt, templtlength, fmm, onsetFunctionLength);
	
	localscorelength =  (templtlength / 2 + onsetFunctionLength - 1) - (templtlength / 2)+1;
	double* localscore =
      slice(localscore_pre, templtlength / 2, templtlength / 2 + onsetFunctionLength - 1);
	free(localscore_pre);

	// backlink(time) is index of best preceding time for this point
	// cumscore(time) is total cumulated score to this point
  int* backlink = (int*)malloc(localscorelength * sizeof(int));
  memset(backlink, 0, localscorelength * sizeof(int)); //necc?
	double* cumscore = (double*)malloc(localscorelength * sizeof(double));
  memset(cumscore, 0, localscorelength * sizeof(int)); //necc?

	// search range for previous beat
	int prangemin = -2 * pd;
	int prangemax = -pd / 2;
	
	// Skewed window
	// txwt = exp(-0.5*((tightness*log(prange/-pd)).^2));
	double tightness = 6.0;
	int txwtlength = (1 + (prangemax - prangemin));
	temp_range = range(prangemin, prangemax, 1);
	temp_times = bdp_times_array_scalar( temp_range ,txwtlength, -1.0 / pd);
	double* temp_log = bdp_log(temp_times, txwtlength);
	temp_times2 = bdp_times_array_scalar(temp_log, txwtlength, tightness);
	temp_power = power(temp_times2, txwtlength, 2.0);
	double* temp_times3 = bdp_times_array_scalar(temp_power, txwtlength, -0.5);
	double* txwt = bdp_exp(temp_times3, txwtlength);
	free(temp_range);
	free(temp_times);
	free(temp_log);
	free(temp_times2);
	free(temp_power);
	free(temp_times3);

	int starting = 1;
	double maxlocalscore = bdp_max(localscore, localscorelength);
	double alpha = 0.9;

	for (int i = 0; i < localscorelength; i++) {
		double* scorecands;// = new double[txwt.length];allocate only in if clause!
		int scorecandslength = txwtlength;//same in both clauses
		
		// Are we reaching back before time zero?
		if (i + prangemin < 0) {
			int valpts = 0;
			scorecands=(double*)malloc(scorecandslength * sizeof(double));
      memset(scorecands, 0, scorecandslength * sizeof(double));

			if (i + prangemax >= 0) {
				valpts = i + prangemax + 1;
        //taken from first slice() below;
				int valvalslength = ((txwtlength - 1)-(txwtlength - valpts)+1);
				double* temp_slice = slice(txwt, txwtlength - valpts, txwtlength - 1);
				double* temp_slice2 = slice(cumscore, i + prangemax- valpts + 1, i + prangemax);
				double* valvals = bdp_times_array_array(temp_slice, valvalslength, temp_slice2);
				free(temp_slice);
				free(temp_slice2);
        for (int j = 0; j < valpts; j++) {
					scorecands[scorecandslength - valpts + j] = valvals[j];
        }
				free(valvals);
			}
		} else {
			// Search over all possible predecessors and apply transition weighting
			double* temp_slice = slice(cumscore, i + prangemin, i + prangemax);
			scorecands = bdp_times_array_array(txwt, txwtlength, temp_slice);
			free(temp_slice);
		}

		// Find best predecessor beat
		double vv = bdp_max(scorecands, scorecandslength);
		int xx = argmax(scorecands, scorecandslength);
		// Add on local score
		cumscore[i] = alpha * vv + (1 - alpha) * localscore[i];
		
		// special case to catch first onset
		if (starting == 1 && localscore[i] < 0.01 * maxlocalscore) {
			backlink[i] = -1;
		} else {
			backlink[i] = i + prangemin + xx;
			// prevent it from resetting, even through a stretch of silence
			starting = 0;
		}
		free(scorecands);
	}
	
	free(txwt);

	long cumscorelength = localscorelength;

	// Backtrace
	// Cumulated score is stabilized to lie in constant range, so just look for one near the end that
  // has a reasonable score
	
  int* templocalmax = localmax(cumscore, cumscorelength);
  subsreflength = 0;
	temp_subref = subsref(cumscore, cumscorelength,templocalmax , cumscorelength, &subsreflength);
  double medscore = median(temp_subref, subsreflength);//tempsublength);
	free(templocalmax);
	free(temp_subref);
	long bestendx = 0;
	long jj = cumscorelength - 2;
	while (jj > 0 && bestendx == 0) {
		if (cumscore[jj] > cumscore[jj - 1] &&
			cumscore[jj] >= cumscore[jj + 1] &&
			cumscore[jj] > 0.5 * medscore) {
			bestendx = jj;
		}
		--jj;
	}
  if (verbose) {
    printf("\n bestendx = %ld medscore = %f", bestendx, medscore);
  }

	int nbeats = 0;
  long* tmplinks = (long*)malloc(cumscorelength * sizeof(long));
	tmplinks[0] = bestendx;
	int bb;
	while ((bb = backlink[tmplinks[nbeats]]) > 0) {
		++nbeats;
		tmplinks[nbeats] = bb;
	}
	
	++nbeats;
	*outputLength = nbeats;
	long* beatData = (long*)malloc(nbeats * sizeof(long));
	
	for (int j = 0; j < nbeats; ++j) {
		beatData[j] = tmplinks[nbeats - 1 - j]; // in frame indexes
	}
	
	// free everything except |beatData|.
	free(fmm);
	free(xfmm);
	free(templt);
	free(localscore);
	free(xpks);
	free(backlink);
	free(cumscore);
  free(tmplinks);

  // end of checkOnsets. Take |data| and do any additional massaging of the data to see if it can
  // run within min/max BPMs.

  if (!beatData || nbeats < 2 || subdivision == 0) {
    printf("\nCould not massage beat data");
    *outputLength = 0;
    return NULL;
  }

  float newBeatToOldBeatRatio = 1; // 2 for twice as many new, .5 for half
  int takeOneOfXBeats = 1; // Used to derive which subdivision to start counting on.
  // Average number of samples between beats. Note that we divide by length-1 since the last beat
  // is not "complete" in terms of having a duration.
  float avgSampleInterval =
      (float)(beatData[nbeats - 1] - beatData[0]) / (nbeats - 1) * _hopSize;
  float avgBPM = (float)60 * _sampleRate / avgSampleInterval;

  // Fold the data to be within the min and max bpms. It is assumed that maxBPM = 2*minBPM, so that
  // duplet feels can always be folded somewhere within that range. Keep folding until in the
  // desired range, but catch edge cases where triplet/quadruplet feels cannot be folded.
  // TODO: no, recursion doesn't work now with the ratios, i.e. shouldn't keep tripling.
  // Assume one fold is enough.
  //while(avgBPM < minBPM || avgBPM > maxBPM) {
  if (avgBPM < minBPM) { // if slow, smaller than min bpm ,
    if (subdivision == 3 && avgBPM * subdivision > maxBPM) { // slow triplet edge case
      if (verbose) {
        printf("\nedge case where avgBPM %f, subdivision %d", avgBPM, subdivision);
      }
      subdivision = 1;// leave as is outside range
    }
    if (subdivision == 4 && avgBPM * subdivision > maxBPM) { // fold quad into duple
      subdivision = 2;
    }

    avgBPM *= subdivision;
    newBeatToOldBeatRatio *= subdivision;
  } else if (avgBPM > maxBPM) { // if fast, greater than max bpm
    // assumably, neither a 2x nor 3x period was not found. All we can do is halve the bpm.
    avgBPM *= .5;
    newBeatToOldBeatRatio *= .5;
    takeOneOfXBeats *= 2;
  }

  if (fabs(newBeatToOldBeatRatio - 1.0f) < FLT_EPSILON) {
    // No change to original data, pass back original values.
    *outputLength = nbeats;
    return beatData;
  }

  long* massgedBeatData = NULL;
  // Halve (or third/quarter/etc) the number of beats.
  if (newBeatToOldBeatRatio < 1.) {
    // Generate small array containing the sum onset values of that subdivision.
    double* beatModVal = (double*)malloc(takeOneOfXBeats * sizeof(double));
    memset(beatModVal, 0, takeOneOfXBeats * sizeof(double));

    // Populate an array of the onset value at each beat.
    float *onsetVals = (float*)malloc(nbeats * sizeof(float));
    for(long i = 0; i < nbeats; i++) {
      onsetVals[i] = onsetFunction[beatData[i]]; //onset function value of the frame of beat i
    }

    // Compute sum of the onset energy of each subdivision.
    for (int i = 0; i < nbeats; i++) {
      beatModVal[i%takeOneOfXBeats] += onsetVals[i];
    }

    free(onsetVals);

    // Find the strongest subdivision index.
    int startBeat = -1;
    double bestBeatVal = -1;
    for (int i = 0; i < takeOneOfXBeats; i++) {
      if (beatModVal[i] > bestBeatVal) {
        bestBeatVal = beatModVal[i];
        startBeat = i;
      }
    }

    free(beatModVal);

    int length = (int)((nbeats - startBeat) * newBeatToOldBeatRatio);
    massgedBeatData = (long*)malloc(length * sizeof(long));
    *outputLength = length;

    for (int i = 0; i < length; i++) {
      massgedBeatData[i] = beatData[i * takeOneOfXBeats + startBeat];
    }
  } else if(newBeatToOldBeatRatio > 1.) { //double or more============
    int stepper = (int)newBeatToOldBeatRatio; //should be 2 or 4

    //so that last original beat is last beat, not adding 1 or three to end...
    int length = (int)(nbeats * newBeatToOldBeatRatio) - (stepper - 1);

    massgedBeatData = (long*)malloc(length * sizeof(long));
    *outputLength = length;

    for(int i = 0, j = 0; i < length; i += stepper, j++) {
      massgedBeatData[i] = beatData[j];
      if(j < nbeats - 1) { //if not last beat
        for(int k = 1; k < stepper; k++) {
          //interpolate
          massgedBeatData[i + k] =
              beatData[j] + (beatData[j + 1] - beatData[j]) * ((float)k / stepper);
        }
      }
    }
  }
  free(beatData);
  return massgedBeatData;
}

// more helper functions

double lin2mel(double fq) {
	return 1127.0 * log(1.0 + (fq / 700.0));
}

double mel2lin(double mel) {
	return 700.0 * (exp(mel / 1127.0) - 1.0);
}

void makeWindow(int size) {
	// Make a Blackman window.
	// w(n)=0.42-0.5cos{(2*PI*n)/(N-1)}+0.08cos{(4*PI*n)/(N-1)};
	_window = (float*)malloc(size * sizeof(float));
  for (int i = 0; i < size; i++) {
		_window[i] = 0.42 - 0.5 * cos(2 * M_PI * i / (size - 1)) + 0.08 * cos(4 * M_PI * i / (size - 1));
  }
}

// Initializer
void AccelInit(int analysisFrameSize, int sampleRate, int hopSize, int melSpecSize) {
  _analysisFrameSize = analysisFrameSize;
  _sampleRate = sampleRate;
  _hopSize = hopSize;
  _melSpecSize = melSpecSize;
  int halfFrameSize = _analysisFrameSize / 2;

  vDSP_Length log2n = (vDSP_Length)log2f(_analysisFrameSize);
  _fftSetup = vDSP_create_fftsetup(log2n, FFT_RADIX2);

  _complexComponents.realp = (float*)malloc(halfFrameSize * sizeof(float));
  _complexComponents.imagp = (float*)malloc(halfFrameSize * sizeof(float));

  _floatShuttle = (float*)malloc(_analysisFrameSize*sizeof(float));

	_melSpec=(float*)malloc((_melSpecSize) * sizeof(float));
	_melSpecD=(double*)malloc((_melSpecSize) * sizeof(double));

  makeWindow(_analysisFrameSize);
	//slower stuff
	_melOfLin = (double*)malloc((halfFrameSize + 1) * sizeof(double));
	_melCenter = (double*)malloc((_melSpecSize + 2) * sizeof(double));
	_melWidth = (double*)malloc((_melSpecSize + 2) * sizeof(double));

	double hzPerBin = (double)sampleRate / 2 / (halfFrameSize + 1);
	double tempMultiplier;
	double melMin, melMax;
	
	//init melCenter
	melMin = lin2mel(0);
	melMax = lin2mel(8000); //assumes non-degenerate sample rate (of > 16K)

  tempMultiplier = (melMax - melMin) / (_melSpecSize + 1);
	for (int i = 0; i < _melSpecSize + 2; i++) {
		_melCenter[i] = melMin + i * tempMultiplier;
	}

	//init MelWidth
	double linbinwidth;
	for (int i = 0; i < _melSpecSize + 1; i++) {
		_melWidth[i] = _melCenter[i + 1] - _melCenter[i];
		linbinwidth = (mel2lin(_melCenter[i + 1]) - mel2lin(_melCenter[i]))/ hzPerBin;
		if (linbinwidth < 1) {
			_melWidth[i] = lin2mel(mel2lin(_melCenter[i]) + hzPerBin)- _melCenter[i];
		}
	}

	//init melOfLin
	tempMultiplier = ((double)sampleRate / (2 * (halfFrameSize + 1)));//same as hzPerBin, right?
	for (int i = 0; i < (halfFrameSize+1); i++) {
		_melOfLin[i] = lin2mel(i*tempMultiplier);
	}
	//end of slower stuff
	
	//Lookup table for weights hard coded 512 bins, 40 melspecdim.
	double tempWeights[100];

  _weightArrayCount = (int *)malloc(_melSpecSize * sizeof(int));
  _weightArrayStartBin = (int *)malloc(_melSpecSize * sizeof(int));
  _ptrArrayToWeightArrays = (double **)malloc(_melSpecSize *sizeof(double*));

	for (int bin = 0; bin < _melSpecSize; bin++) {
		_weightArrayStartBin[bin] = -1;
		int tempWeightIndex = 0;
		for (int i = 0; i < halfFrameSize; i++) {
			double weight = 1.0 - (fabs(_melOfLin[i] - _melCenter[bin]) / _melWidth[bin]);
			if (weight > 0) {
				if (_weightArrayStartBin[bin] == -1) {
					_weightArrayStartBin[bin] = i;
				}
				tempWeights[tempWeightIndex]=weight;
				tempWeightIndex++;
			}
		}
		//for this melspec bin, create little array of weights, and remember the size
		_weightArrayCount[bin] = tempWeightIndex;
		_ptrArrayToWeightArrays[bin] = (double*)malloc(tempWeightIndex * sizeof(double));
		for (int i=0; i < tempWeightIndex; i++) {
			_ptrArrayToWeightArrays[bin][i] = tempWeights[i];
		}
	}
  _rms = rmsTarget;
}

void AccelDealloc() {
  vDSP_destroy_fftsetup(_fftSetup);
  free(_complexComponents.imagp);
  free(_complexComponents.realp);
  free(_floatShuttle);
  free(_window);
  free(_melSpec);
  free(_melSpecD);
  free(_melOfLin);
  free(_melCenter);
  free(_melWidth);
  free(_weightArrayCount);
  free(_weightArrayStartBin);

  for (int i = 0; i < _melSpecSize; i++) {
    free(_ptrArrayToWeightArrays[i]);
  }
  free(_ptrArrayToWeightArrays);
}

float* AccelRender(short* ioData) {
  int halfAnalysisFrameSize = _analysisFrameSize >> 1;

	// Convert short input values to array of float.
	vDSP_vflt16(ioData, 1, _floatShuttle, 1, _analysisFrameSize);

  //scalar for division to normalize values to 0-1.
  //float *shortToFloatDivisor = (float*)malloc(1 * sizeof(float));
  //shortToFloatDivisor[0] = (float)SHRT_MAX;
  float shortToFloatDivisor = (float)SHRT_MAX;
	vDSP_vsdiv(_floatShuttle, 1, &shortToFloatDivisor , _floatShuttle, 1, _analysisFrameSize);
  //free(shortToFloatDivisor);
	
	// Normalize audio's RMS using a moving average estimate of it.
	// Calculate current frame's RMS:
	float currentFrameRMS = 0;
	vDSP_rmsqv(_floatShuttle, 1 ,&currentFrameRMS, _analysisFrameSize);//rms func does sqrt( 1/n * sum floatshuttle[i]^2
	
	// Update moving average
	_rms = rmsAlpha * currentFrameRMS + (1 - rmsAlpha) * _rms;

	// Mult by window vector
  vDSP_vmul (_floatShuttle,1,_window,1,_floatShuttle,1,_analysisFrameSize );

	//Mult by rmsTarget/rms scalar.
	float rmsMultiplier = rmsTarget/_rms; //DEI safety., clip at 1.
  //printf("\nrms: %.2f", rmsMultiplier);
	vDSP_vsmul (_floatShuttle, 1, &rmsMultiplier, _floatShuttle, 1, _analysisFrameSize);
	
	//Convert real input to even-odd components.
  vDSP_ctoz((DSPComplex *)_floatShuttle, 2, &_complexComponents, 1, halfAnalysisFrameSize);

  //Setup and run FFT.
  float ln = log2f(_analysisFrameSize); //DEI check number
  vDSP_fft_zrip(_fftSetup, &_complexComponents, 1, ln, FFT_FORWARD);

  // Get magnitudes and store in _complexComponents.realp
  // Magnitudes are real*real+imag*imag, but _not_ square rooted!
  vDSP_zvmags(&_complexComponents, 1, _complexComponents.realp, 1, halfAnalysisFrameSize);

  // Normalize.
  float scale = 1.f/(halfAnalysisFrameSize);
  vDSP_vsmul(_complexComponents.realp, 1, &scale, _complexComponents.realp, 1, halfAnalysisFrameSize);

  // Convert the FFT magnitudes into a melodic spectrum.
	for (int bin = 0; bin < _melSpecSize; bin++) {
		_melSpecD[bin] = 0;
		for (int i=0;i<_weightArrayCount[bin];i++) {
			_melSpecD[bin]+=_ptrArrayToWeightArrays[bin][i] * _complexComponents.realp[ _weightArrayStartBin[bin]+i ];
		}
	}

  // Set to decibel representation
	double one = 1;
	vDSP_vdbconD(_melSpecD, 1, &one, _melSpecD, 1, _melSpecSize, 0);
  for (int bin = 0; bin < _melSpecSize; bin++) {
    if (_melSpecD[bin]<melSpecSilenceThresh) {
      _melSpecD[bin]=melSpecSilenceThresh;
    }
  }
	
  for (int bin = 0; bin < _melSpecSize; bin++) {
    _melSpec[bin]=(float)_melSpecD[bin];
  }

	return _melSpec;
}
