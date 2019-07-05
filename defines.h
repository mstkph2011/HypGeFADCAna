//global defines for the analysis, add this header file to any file of the analysis

#ifndef DEFINES
#define DEFINES

//#define TRACE_LENGTH 10
//#define TRACE_LENGTH 2048
#define TRACE_LENGTH 4096
//#define TRACE_LENGTH 16384
//#define TRACE_LENGTH 32768
#define FADC_CHAN 1
// second one is actuaclly NaI!!!!!! taken care of in AnlProc files, changed to 1 there. Change for Jülich 2015 (steinen, 06.05.2015)
// FADC_CHAN max 8 !!	

// this factor is calculated by the product of the Sampling Rate and the wanted x axis scale of your trace histograms. e.g. 100 MSa/s * 10⁻6 s --> 0.01 , axis in \mu s
#define TIME_RESOLUTION_FACTOR 0.01
#endif //DEFINES
