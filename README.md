## Real-time SDR for mono/stereo FM and RDS

The main objective of the project was to navigate a complex speciﬁcation and understand the challenges that must be addressed for a real-time implementation of a computing system that operates in a form factor-constrained environment. 

## Repository Structure Breakdown
 
Below is a quick breakdown of what each folder is used for.
  - `src`
    - This folder is where all the C++ Source Code is written. Running the radio executable with either captured Raw Data or piping in output from an SDR was used         to test throughout development.
  - `model`
    - This folder is where all of our unoptimized python models reside. Some numpy libraries were used which were later optimized to run in real time with the SDR
    -  Signal Flow Graphs were tested here first and then further implemented and optmized in C++
  - `doc`
    - Contains Report evaluating performance and overall implementation of signal flow graph for the mono/stereo FM standard
  - `data`
    - Empty..Used to contain RAW files used during development but due to space limitations was removed. Can be either downloaded online or be captured using an SDR
    
    
