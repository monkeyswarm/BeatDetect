# BeatDetect
Beat-based audio analysis and playback.

BeatDetect is a tool for audio beat detection and tempo-based playback. It is a (currently OSX only) command-line program, which analyzes an audio file for beat positions, and (optionally) plays the audio back at a user-specified BPM. The command-line wrapper and audio playback are OSX-specific (in Objective C), but the core analysis routine is in C. Future work will make the whole project more cross-platform.

## Example usage:
This tool uses OSX-style flags as the command line parameters.
```
// Analyze and play foo.mp3 at 120 BPM. Add clicks on the beats.
$ ./BeatDetect -file foo.mp3 -BPM 120 -shouldClick YES

// Analyze foo.mp3 and write beat data to a text file.
$ ./BeatDetect -file foo.mp3 -shouldPlay NO -outputFile output.txt
```

## About the analysis and playback:
The analysis algorithm is based on Dan Ellis' [dynamic programming algorithm](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.94.8773&rep=rep1&type=pdf), but with my own modifications that nudge the resulting BPM, when appropriate, to within a desired range (currently 100 to 200 BPM). The algorithm expects an overall constant tempo, though will successfully track individual beats as the tempo varies. 

The analysis:
1) Takes decompressed audio data, runs an FFT, generating a mel spectrum at each frame.
2) Uses at the differences between spectra as an onset function.
3) Uses autocorrelation to detect overall periodicity in the onset function.
4) Chooses an overall periodicity (and relevant subdivisions, e.g. duple, triplet, etc) based on the autocorrelation peaks and the desired output BPM range.
5) Uses dynamic programming to maximize individual beat energy and adherence to an overall tempo.
6) If result is still outside desired BPM range, see if it can be "folded" based on its subdivisions. E.g. a 60 BPM result, which had a strong periodicity at twice the speed (but not at three times the speed) could be recast as 120 BPM.
7) Returns the beat map of positions, in seconds.

The playback algorithm is effectively a time-domain time-stretcher, but with window length based on beat position and playback rate.
Based on a user-selected BPM, a phasor runs at that frequency. It snaps playback position to the relevant beat and subdivision 4 times per beat (i.e. on 16th notes). A small crossfade resolves discontinuities between windows. With the windows synchronized the the tempo, aural artifacts are far less noticable.

## Originally from...
This is the engine of the (currently removed from sale) iOS app [MiniMash](https://vimeo.com/25455527)

## License
Released under GNU GPL v3. See LICENSE.txt.

### Contact
dan -at- danieliglesia.com
