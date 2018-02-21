#include "Util.h"
#include <map>
#include <string>
#include <cfloat>
#include <cstdlib>
#pragma once
using namespace std;


struct TimeInfo{
    std::string name;
    bool started;
    int64 run_time;
    float total_time_seconds;
    float max_time_seconds;
    float min_time_seconds;
    float average_call_time;
    float self_timer;
    int calls;
};


class Timer
{
    public:
        static Timer& getInstance()
        {
            static Timer    instance; // Guaranteed to be destroyed.
                                  // Instantiated on first use.
            return instance;
        }
        
        void StartTimer(std::string name);
        void EndTimer(std::string name);
        float GetTotalTimeInSeconds(std::string name);
        float GetAvgTimeInSeconds(std::string name);
        float GetMaxTimeInSeconds(std::string name);
        float GetMinTimeInSeconds(std::string name);
        int GetCallCount(std::string name);
        void PrintTimerInfo();
        
        
        
        
    private:
        Timer() {
        };                   // Constructor? (the {} brackets) are needed here.
        
        std::map<std::string,TimeInfo> time_data;
        
        bool timer_exists(const std::string & name);
        
    void initialize_timer(std::string name);
    //Un-implemented so you cant mess with this singleton :)
        Timer(Timer const&);              // Don't Implement
        void operator=(Timer const&); // Don't implement
};

