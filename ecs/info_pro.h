#ifndef __INFO__PRO__
#define __INFO__PRO__

#include<string>
#include<vector>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <numeric>
#include "lib_io.h"
#include <iostream>
#include <sstream>
#include <ctime>
#include <iomanip>
#include "LMatrix.h"

const int oneDayLong = 86400;
using namespace std;

extern int predict_daySpan, history_daySpan; //需要预测的时间跨度和总的历史时间跨度
extern vector<int> vflavors_pridict_nums;  //整理好的预测值标签
extern string dim, predict_begin_time, predict_end_time, history_begin_time, history_end_time;
extern vector< vector<int> > sequences;
extern time_t predict_begin_time_t, predict_end_time_t, history_begin_time_t, history_end_time_t;

struct PhyInfo {
	int phycpu;
	int phymem;
	int phyhard;
} ;
extern PhyInfo phyinfo;

struct FlavorsInfo {
	vector<string> vflavors;
	vector<int> vcpus, vmems;
	vector<double> bizong;
} ;
extern FlavorsInfo flavorsinfo;


//int oneDayLong = 24 * 3600;

void get_info(char * info[MAX_INFO_NUM]);

void get_data(char * data[MAX_DATA_NUM], int data_num);

time_t str2time(string str, string format);

void simple_predict();
void lwlr_predict();
void aclwlr_predict();
void gray_lwlr_predict();
void cig_lwlr_predict();
void ThreeIndex_predict();




#endif