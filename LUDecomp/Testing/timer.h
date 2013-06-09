#ifndef TIMER_H
#define TIMER_H 1

#include <sys/time.h>

class Timer{

private:
    timeval _start;
    timeval _stop;
    
public:
    
    void
    start(){
        gettimeofday(&_start, 0);
    }
    
    void
    stop(){
        gettimeofday(&_stop, 0);
    }
    
    double
    elapsed(){
        long long int usecs_start = _start.tv_sec * 1000000. + _start.tv_usec;
        long long int usecs_stop =  _stop.tv_sec * 1000000. + _stop.tv_usec;
        return double(usecs_stop - usecs_start) / 1000000.;
    }
};

#endif // TIMER_H