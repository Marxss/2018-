#include "predict.h"
#include <stdio.h>
#include "allocate.h"



//��Ҫ��ɵĹ��������
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
	finish = clock(); // �˻���̽���
	double duration = ((double)(finish - start)) / CLOCKS_PER_SEC;
	cout << "�����㷨��ʱ��" << duration << endl;
	// ��Ҫ���������
	//char * result_file = (char *)"17\n\n0 8 0 20";

	// ֱ�ӵ�������ļ��ķ��������ָ���ļ���(ps��ע���ʽ����ȷ�ԣ�����н⣬��һ��ֻ��һ�����ݣ��ڶ���Ϊ�գ������п�ʼ���Ǿ�������ݣ�����֮����һ���ո�ָ���)
	//write_result(result_file, filename);
}
