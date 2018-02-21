#include "Timer.h"

#pragma once
/*---------------Public METHODS--------------------------*/  

int sort_time_func (const void * a, const void * b)
{
    TimeInfo* tima,*timb;
    tima = (TimeInfo*)a;
    timb = (TimeInfo*)b;
//    printf("nu a = %11.6f nu b = %11.6f\n",nu_a,nu_b);
  if (tima->total_time_seconds >  timb->total_time_seconds ) return -1;
  if ( tima->total_time_seconds == timb->total_time_seconds) return 0;
     return 1;
};


void Timer::StartTimer(std::string name){
    initialize_timer(name);
    
    if(time_data[name].started){
        fprintf(stderr,"[Timer::StartTimer]: Timer of name %s has already started!!\n",name.c_str());
        //exit(0);
    }else{
        time_data[name].started = true;
        time_data[name].run_time = GetTimeMs64();
        time_data[name].calls++;
    }

};
void Timer::EndTimer(std::string name){

    if( !timer_exists(name)){
        fprintf(stderr,"[Timer::EndTimer]: Timer of name %s does not exist!!\n",name.c_str());
        return;
        //exit(0);    
    }
    TimeInfo temp = time_data[name];
    if(!temp.started){
        fprintf(stderr,"[Timer::EndTimer]: Timer of name %s has not been started!!\n",name.c_str());
        //exit(0);
    }else{
        
        float time_taken_seconds =float(GetTimeMs64()- temp.run_time)/1000.0f;
        temp.started = false;
        temp.total_time_seconds += time_taken_seconds;
        temp.max_time_seconds = max(temp.max_time_seconds,time_taken_seconds);
        temp.min_time_seconds = min(temp.min_time_seconds,time_taken_seconds);
        temp.average_call_time =  temp.total_time_seconds/(float)temp.calls;
        time_data[name] = temp;
        for (std::map<std::string,TimeInfo>::iterator it=time_data.begin(); it!=time_data.end(); ++it){
            if(it->second.started)
                it->second.self_timer -= time_taken_seconds;
        }
    }
};


        
float Timer::GetTotalTimeInSeconds(std::string name){
    if(!timer_exists(name)){
        fprintf(stderr,"[Timer::GetTotalTimeInSeconds]: Timer of name %s does not exist!!\n",name.c_str());
        exit(0);    
    }
    
    if(time_data[name].started){
        float saved_time = time_data[name].total_time_seconds;
        saved_time += (GetTimeMs64() - time_data[name].run_time)/1000.0f;
        return saved_time;
    }
    
    return time_data[name].total_time_seconds;
};
float Timer::GetAvgTimeInSeconds(std::string name){
    if(!timer_exists(name)){
        fprintf(stderr,"[Timer::GetAvgTimeInSeconds]: Timer of name %s does not exist!!\n",name.c_str());
        exit(0);    
    }
    return time_data[name].average_call_time;
};
float Timer::GetMaxTimeInSeconds(std::string name){
    if(!timer_exists(name)){
        fprintf(stderr,"[Timer::GetMaxTimeInSeconds]: Timer of name %s does not exist!!\n",name.c_str());
        exit(0);    
    }
    return time_data[name].max_time_seconds;
};
float Timer::GetMinTimeInSeconds(std::string name){    
    if(!timer_exists(name)){
        fprintf(stderr,"[Timer::GetMinTimeInSeconds]: Timer of name %s does not exist!!\n",name.c_str());
        exit(0);    
    }
    return time_data[name].min_time_seconds;

};

int Timer::GetCallCount(std::string name){
    if(!timer_exists(name)){
        fprintf(stderr,"[Timer::GetAvgTimeInSeconds]: Timer of name %s does not exist!!\n",name.c_str());
        exit(0);    
    }
    return time_data[name].average_call_time;
};

        
void Timer::PrintTimerInfo(){

    if(time_data.empty()){
        printf("[Timer::PrintTimerInfo]: No timer data collected!!\n");
        return;
    }
    
    //Otherwise we get the count
    int timer_size = time_data.size();
    //Allocate array of that size
    TimeInfo* temp = new TimeInfo[timer_size];
    int i = 0;
      for (std::map<std::string,TimeInfo>::iterator it=time_data.begin(); it!=time_data.end(); ++it){
            
            if(it->second.started){
            EndTimer(it->first);
            //float time_taken_seconds =float(GetTimeMs64()- temp[i].run_time)/1000.0f;
            //temp[i].started = false;
            //temp[i].total_time_seconds += time_taken_seconds;
            //temp[i].max_time_seconds = max(temp[i].max_time_seconds,time_taken_seconds);
            //temp[i].min_time_seconds = min(temp[i].min_time_seconds,time_taken_seconds);
            //temp[i].average_call_time =  temp[i].total_time_seconds/(float)temp[i].calls;                
            }
        //EndTimer(it->first);
            temp[i] = it->second;
            
            i++;
        }
        
        qsort(temp,timer_size,sizeof(TimeInfo),sort_time_func);
        //Print out timer information
        printf("--------------------------------------------------------------------------------------------\n");
        printf("--                                   Timer Data                                           --\n");
        printf("--------------------------------------------------------------------------------------------\n");
        printf("--                      Time data sorted by longest total time                            --\n");
        printf("--------------------------------------------------------------------------------------------\n");
        printf("--------------------------------------------------------------------------------------------\n");
        printf("Name                    Calls      Total time         Avg        Max        Min     Self    \n");   
         printf("--------------------------------------------------------------------------------------------\n");       
        for(int i = 0; i < timer_size; i++)
            printf("%18.18s   %8d   %12.6fs     %6.4fs     %6.4fs     %6.4fs   %6.4fs\n",temp[i].name.c_str(),temp[i].calls
                                                           ,temp[i].total_time_seconds
                                                           ,temp[i].average_call_time,temp[i].max_time_seconds
                                                           ,temp[i].min_time_seconds,temp[i].self_timer+temp[i].total_time_seconds);
                

        printf("---------------------------------------------------------------------------\n");
        printf("---------------------------------------------------------------------------\n");
        
        delete[] temp;




};
        
        
        
/*---------------PRIVATE METHODS--------------------------*/        


bool Timer::timer_exists(const std::string & name){return time_data.count(name);}     


void Timer::initialize_timer(const std::string name){

    if( !timer_exists(name) ){
        time_data[name].name = name;
        time_data[name].started = false;
        time_data[name].run_time = 0;
        time_data[name].total_time_seconds = 0.0f;
        time_data[name].max_time_seconds = 0.0f;
        time_data[name].min_time_seconds = FLT_MAX;
        time_data[name].average_call_time = 0.0f;
        time_data[name].calls = 0;
        time_data[name].self_timer = 0.0;

    }

};

