#ifndef TIMER_H
#define TIMER_H 1

#include <sys/time.h>

class Timer{

private:
    timeval _start;
    timeval _stop;
    
    bool _is_start;
    
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
        return ((_stop.tv_sec - _start.tv_sec)*1000000u + (_stop.tv_usec-_start.tv_usec))/1000000.0;
    }
};

#endif // TIMER_H