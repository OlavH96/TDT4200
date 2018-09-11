#include <chrono>
using namespace std::chrono;

class Timer {

public:
  long long int getIterations(){
    return iterations;
  }
  long long int getTotalTime(){
    return totalTime;
  }
  void iterate(){
    timer_end = high_resolution_clock::now();
    totalTime += getTimeElapsedMilliseconds();
    iterations++;
    resetStart();
  }
  long long int getAverageTime(){
    return totalTime / iterations;
  }
  void resetStart(){
    if(running){
      timer_start = high_resolution_clock::now();
    }
  }
  void start(){
    timer_start = high_resolution_clock::now();
    running = true;
  }
  void stop(){
    timer_end = high_resolution_clock::now();
    running = false;
  }
  long long int getTimeElapsedMilliseconds(){
    if(running){
      timer_end = high_resolution_clock::now();
      return duration_cast<nanoseconds>(timer_end - timer_start).count();
    }
    else{
      return duration_cast<nanoseconds>(timer_end - timer_start).count();
    }
  }
  bool isRunning(){
    return running;
  }
private:
  high_resolution_clock::time_point timer_start;
  high_resolution_clock::time_point timer_end;
  bool running = false;
  long long int totalTime = 0;
  long long int iterations = 0;
};
