#include "predict.h"
#include <stdio.h>
#include "allocate.h"



//你要完成的功能总入口
void predict_server(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, char * filename)
{
	time_t start, finish;
	start = clock();
	get_info(info);

	get_data(data, data_num);
	//simple_predict();
	lwlr_predict();
	//cig_lwlr_predict();
	//aclwlr_predict();
	//gray_lwlr_predict();
	//ThreeIndex_predict();

	allocate_simple(flavorsinfo.vflavors, vflavors_pridict_nums, dim, filename);
	finish = clock(); // 退火过程结束
	double duration = ((double)(finish - start)) / CLOCKS_PER_SEC;
	cout << "放置算法用时：" << duration << endl;
	// 需要输出的内容
	//char * result_file = (char *)"17\n\n0 8 0 20";

	// 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
	//write_result(result_file, filename);
}
