#include "allocate.h"
#include <algorithm>
struct PhyInfo {
	int phycpu;
	int phymem;
	int phyhard;
};
extern PhyInfo phyinfo;

struct FlavorsInfo {
	vector<string> vflavors;
	vector<int> vcpus, vmems;
};
extern FlavorsInfo flavorsinfo;

struct deploy {
	int phyid;
	vector<string> vflavors;
	vector<int> nums;  //能放置该类型虚拟机个数
	int cpuLeft;
	int memLeft;
};

void allocate_simple(vector<string> vflavors, vector<int> vflavors_pridict_nums, string dim, char* filename)
{

	//三类排序法
	vector<int> vflavors_pridict_nums_copy = vflavors_pridict_nums;
	vector<struct deploy> deploys;
	deploy new_deploy;
	new_deploy.phyid = deploys.size() + 1;
	vector<int> nums(flavorsinfo.vflavors.size(), 0);
	new_deploy.nums = nums;
	new_deploy.cpuLeft = phyinfo.phycpu;
	new_deploy.memLeft = phyinfo.phymem;

	FlavorsInfo small_vm;
	FlavorsInfo middle_vm;
	FlavorsInfo large_vm;

	for (int i = flavorsinfo.vflavors.size() - 1; i >= 0; i--)
	{

		int mems_cpus = flavorsinfo.vmems[i] / flavorsinfo.vcpus[i];
		if (mems_cpus == 4096)
		{
			while (vflavors_pridict_nums[i])
			{
				large_vm.vflavors.push_back(flavorsinfo.vflavors[i]);
				large_vm.vcpus.push_back(flavorsinfo.vcpus[i]);
				large_vm.vmems.push_back(flavorsinfo.vmems[i]);

				vflavors_pridict_nums[i]--;
			}
		}
		if (mems_cpus == 2048)
		{
			while (vflavors_pridict_nums[i])
			{
				middle_vm.vflavors.push_back(flavorsinfo.vflavors[i]);
				middle_vm.vcpus.push_back(flavorsinfo.vcpus[i]);
				middle_vm.vmems.push_back(flavorsinfo.vmems[i]);
				vflavors_pridict_nums[i]--;
			}
		}
		if (mems_cpus == 1024)
		{
			while (vflavors_pridict_nums[i])
			{
				small_vm.vflavors.push_back(flavorsinfo.vflavors[i]);
				small_vm.vcpus.push_back(flavorsinfo.vcpus[i]);
				small_vm.vmems.push_back(flavorsinfo.vmems[i]);
				vflavors_pridict_nums[i]--;
			}
		}

	}

	while (large_vm.vflavors.size() || middle_vm.vflavors.size() || small_vm.vflavors.size())
	{
		if (large_vm.vflavors.size() && (large_vm.vcpus[0] <= new_deploy.cpuLeft) && (large_vm.vmems[0] <= new_deploy.memLeft))
		{
			new_deploy.cpuLeft = new_deploy.cpuLeft - large_vm.vcpus[0];
			new_deploy.memLeft = new_deploy.memLeft - large_vm.vmems[0];
			vector<string>::iterator resoult_large = find(new_deploy.vflavors.begin(), new_deploy.vflavors.end(), large_vm.vflavors[0]);
			if (resoult_large == new_deploy.vflavors.end())
			{
				new_deploy.vflavors.push_back(large_vm.vflavors[0]);
				new_deploy.nums[new_deploy.vflavors.size() - 1]++;
			}
			else
			{
				int nPosition_large = distance(new_deploy.vflavors.begin(), resoult_large);
				new_deploy.nums[nPosition_large]++;
			}
			large_vm.vflavors.erase(large_vm.vflavors.begin());
			large_vm.vcpus.erase(large_vm.vcpus.begin());
			large_vm.vmems.erase(large_vm.vmems.begin());
		}
		else
		{
			if (large_vm.vflavors.size() && (large_vm.vcpus[large_vm.vcpus.size() - 1] <= new_deploy.cpuLeft) && (large_vm.vmems[large_vm.vmems.size() - 1] <= new_deploy.memLeft))
			{
				new_deploy.cpuLeft = new_deploy.cpuLeft - large_vm.vcpus[large_vm.vcpus.size() - 1];
				new_deploy.memLeft = new_deploy.memLeft - large_vm.vmems[large_vm.vmems.size() - 1];
				vector<string>::iterator resoult_large = find(new_deploy.vflavors.begin(), new_deploy.vflavors.end(), large_vm.vflavors[large_vm.vflavors.size() - 1]);

				if (resoult_large == new_deploy.vflavors.end())
				{
					new_deploy.vflavors.push_back(large_vm.vflavors[large_vm.vflavors.size() - 1]);
					new_deploy.nums[new_deploy.vflavors.size() - 1]++;
				}
				else
				{
					int nPosition_large = distance(new_deploy.vflavors.begin(), resoult_large);
					new_deploy.nums[nPosition_large]++;
				}
				large_vm.vflavors.pop_back();
				large_vm.vcpus.pop_back();
				large_vm.vmems.pop_back();
			}
			else
			{
				if (middle_vm.vflavors.size() != 0 || small_vm.vflavors.size() != 0)
				{
				}
				else
				{
					if (large_vm.vflavors.size() != 0)
					{
						new_deploy.phyid = deploys.size() + 1;
						deploys.push_back(new_deploy);
						new_deploy.cpuLeft = phyinfo.phycpu;
						new_deploy.memLeft = phyinfo.phymem;
						new_deploy.nums = nums;
						new_deploy.vflavors.clear();
					}
				}
			}
		}
		if (middle_vm.vflavors.size() && (middle_vm.vcpus[0] <= new_deploy.cpuLeft) && (middle_vm.vmems[0] <= new_deploy.memLeft))
		{
			new_deploy.cpuLeft = new_deploy.cpuLeft - middle_vm.vcpus[0];
			new_deploy.memLeft = new_deploy.memLeft - middle_vm.vmems[0];
			vector<string>::iterator resoult_middle = find(new_deploy.vflavors.begin(), new_deploy.vflavors.end(), middle_vm.vflavors[0]);
			if (resoult_middle == new_deploy.vflavors.end())
			{
				new_deploy.vflavors.push_back(middle_vm.vflavors[0]);
				new_deploy.nums[new_deploy.vflavors.size() - 1]++;
			}
			else
			{
				int nPosition_middle = distance(new_deploy.vflavors.begin(), resoult_middle);
				new_deploy.nums[nPosition_middle]++;
			}
			middle_vm.vflavors.erase(middle_vm.vflavors.begin());
			middle_vm.vcpus.erase(middle_vm.vcpus.begin());
			middle_vm.vmems.erase(middle_vm.vmems.begin());
		}
		else
		{
			if (middle_vm.vflavors.size() && (middle_vm.vcpus[middle_vm.vcpus.size() - 1] <= new_deploy.cpuLeft) && (middle_vm.vmems[middle_vm.vmems.size() - 1] <= new_deploy.memLeft))
			{
				new_deploy.cpuLeft = new_deploy.cpuLeft - middle_vm.vcpus[middle_vm.vcpus.size() - 1];
				new_deploy.memLeft = new_deploy.memLeft - middle_vm.vmems[middle_vm.vmems.size() - 1];
				vector<string>::iterator resoult_large = find(new_deploy.vflavors.begin(), new_deploy.vflavors.end(), middle_vm.vflavors[middle_vm.vflavors.size() - 1]);
				if (resoult_large == new_deploy.vflavors.end())
				{
					new_deploy.vflavors.push_back(middle_vm.vflavors[middle_vm.vflavors.size() - 1]);
					new_deploy.nums[new_deploy.vflavors.size() - 1]++;
				}
				else
				{
					int nPosition_large = distance(new_deploy.vflavors.begin(), resoult_large);
					new_deploy.nums[nPosition_large]++;
				}
				middle_vm.vflavors.pop_back();
				middle_vm.vcpus.pop_back();
				middle_vm.vmems.pop_back();
			}
			else
			{
				if (small_vm.vflavors.size() != 0)
				{
				}
				else
				{
					if (middle_vm.vflavors.size() != 0)
					{
						new_deploy.phyid = deploys.size() + 1;
						deploys.push_back(new_deploy);
						new_deploy.cpuLeft = phyinfo.phycpu;
						new_deploy.memLeft = phyinfo.phymem;
						new_deploy.nums = nums;
						new_deploy.vflavors.clear();
					}
				}

			}
		}
		if (small_vm.vflavors.size() && (small_vm.vcpus[0] <= new_deploy.cpuLeft) && (small_vm.vmems[0] <= new_deploy.memLeft))
		{
			new_deploy.cpuLeft = new_deploy.cpuLeft - small_vm.vcpus[0];
			new_deploy.memLeft = new_deploy.memLeft - small_vm.vmems[0];
			vector<string>::iterator resoult_small = find(new_deploy.vflavors.begin(), new_deploy.vflavors.end(), small_vm.vflavors[0]);
			if (resoult_small == new_deploy.vflavors.end())
			{
				new_deploy.vflavors.push_back(small_vm.vflavors[0]);
				new_deploy.nums[new_deploy.vflavors.size() - 1]++;
			}
			else
			{
				int nPosition_small = distance(new_deploy.vflavors.begin(), resoult_small);
				new_deploy.nums[nPosition_small]++;
			}
			small_vm.vflavors.erase(small_vm.vflavors.begin());
			small_vm.vcpus.erase(small_vm.vcpus.begin());
			small_vm.vmems.erase(small_vm.vmems.begin());
		}
		else
		{
			if (small_vm.vflavors.size() && (small_vm.vcpus[small_vm.vcpus.size() - 1] <= new_deploy.cpuLeft) && (small_vm.vmems[small_vm.vmems.size() - 1] <= new_deploy.memLeft))
			{
				new_deploy.cpuLeft = new_deploy.cpuLeft - small_vm.vcpus[small_vm.vcpus.size() - 1];
				new_deploy.memLeft = new_deploy.memLeft - small_vm.vmems[small_vm.vmems.size() - 1];
				vector<string>::iterator resoult_large = find(new_deploy.vflavors.begin(), new_deploy.vflavors.end(), small_vm.vflavors[small_vm.vflavors.size() - 1]);
				if (resoult_large == new_deploy.vflavors.end())
				{
					new_deploy.vflavors.push_back(small_vm.vflavors[small_vm.vflavors.size() - 1]);
					new_deploy.nums[new_deploy.vflavors.size() - 1]++;
				}
				else
				{
					int nPosition_large = distance(new_deploy.vflavors.begin(), resoult_large);
					new_deploy.nums[nPosition_large]++;
				}
				small_vm.vflavors.pop_back();
				small_vm.vcpus.pop_back();
				small_vm.vmems.pop_back();
			}
			else
			{
				if (small_vm.vflavors.size() != 0)
				{
					new_deploy.phyid = deploys.size() + 1;
					deploys.push_back(new_deploy);
					new_deploy.cpuLeft = phyinfo.phycpu;
					new_deploy.memLeft = phyinfo.phymem;
					new_deploy.nums = nums;
					new_deploy.vflavors.clear();
				}
			}
		}
	}

	if (new_deploy.vflavors.size())
	{
		new_deploy.phyid = deploys.size() + 1;
		deploys.push_back(new_deploy);
	}

	/*//将近无穷法
	vector<int> vflavors_pridict_nums_copy = vflavors_pridict_nums;
	vector<int> vflavors_pridict_nums_copy2 = vflavors_pridict_nums;
	vector<int> vflavors_pridict_nums_copy3 = vflavors_pridict_nums;
	vector<int> vflavors_pridict_nums_copy4 = vflavors_pridict_nums;
	vector<struct deploy> deploys;
	vector<struct deploy> deploys2;
	deploy new_deploy;
	new_deploy.phyid = deploys.size() + 1;
	new_deploy.vflavors = flavorsinfo.vflavors;
	vector<int> nums(flavorsinfo.vflavors.size(), 0);
	new_deploy.nums = nums;
	new_deploy.cpuLeft = phyinfo.phycpu;
	new_deploy.memLeft = phyinfo.phymem;
	deploys2.push_back(new_deploy);

	deploy a_new_deploy = deploys2.back();
	int id = -1;

	while (accumulate(vflavors_pridict_nums_copy4.begin(), vflavors_pridict_nums_copy4.end(), 0))
	{
	vflavors_pridict_nums_copy3 = vflavors_pridict_nums_copy4;
	for (int i = 0; i < flavorsinfo.vflavors.size(); i++)
	{
	while (vflavors_pridict_nums_copy3[i])
	{
	vflavors_pridict_nums_copy2 = vflavors_pridict_nums_copy3;
	for (int j = i; j <flavorsinfo.vflavors.size(); j++)
	{
	while ((flavorsinfo.vcpus[j] <= a_new_deploy.cpuLeft) && (flavorsinfo.vmems[j] <= a_new_deploy.memLeft) && (vflavors_pridict_nums_copy2[j]>0))//装满
	{
	a_new_deploy.nums[j]++;
	a_new_deploy.cpuLeft -= flavorsinfo.vcpus[j];
	a_new_deploy.memLeft -= flavorsinfo.vmems[j];
	vflavors_pridict_nums_copy2[j]--;
	}//装满
	}
	if (a_new_deploy.cpuLeft <= deploys2.back().cpuLeft)//如果大于上一个则入栈
	{
	if (a_new_deploy.cpuLeft == deploys2.back().cpuLeft)
	{
	if (a_new_deploy.memLeft < deploys2.back().memLeft)
	{
	deploys2.pop_back();
	deploys2.push_back(a_new_deploy);
	}
	}
	else
	{
	deploys2.pop_back();
	deploys2.push_back(a_new_deploy);
	}
	}//如果大于上一个则入栈
	a_new_deploy.phyid = deploys.size() + 1;
	a_new_deploy.vflavors = flavorsinfo.vflavors;
	vector<int> nums(flavorsinfo.vflavors.size(), 0);
	a_new_deploy.nums = nums;
	a_new_deploy.cpuLeft = phyinfo.phycpu;
	a_new_deploy.memLeft = phyinfo.phymem;
	vflavors_pridict_nums_copy3[i]--;
	}
	}
	deploys.push_back(deploys2.back());
	id++;
	deploys2.push_back(new_deploy);
	for (int a = 0; a < flavorsinfo.vflavors.size(); a++)
	{
	if (vflavors_pridict_nums_copy4[a])
	{
	vflavors_pridict_nums_copy4[a] = vflavors_pridict_nums_copy4[a] - deploys[id].nums[a];
	}
	}
	}
	deploys.back().phyid = deploys.size();
	*/


	//OUTPUT
	string result_file = "";
	//vm nums
	int total_vm = accumulate(vflavors_pridict_nums_copy.begin(), vflavors_pridict_nums_copy.end(), 0);
	char s[20];
	sprintf(s, "%d", total_vm);
	result_file.append(s);
	result_file.append("\n");
	//vm items
	for (int i = 0; i < vflavors.size(); i++)
	{
		char s[20];
		sprintf(s, "%d", vflavors_pridict_nums_copy[i]);
		result_file.append(vflavors[i]).append(" ").append(s).append("\n");
	}
	result_file.append("\n");

	//deploys
	sprintf(s, "%d", deploys.size());
	result_file.append(s).append("\n");
	for (int i = 0; i < deploys.size(); i++)
	{
		char s[20];
		sprintf(s, "%d", deploys[i].phyid);
		result_file.append(s).append(" ");
		for (int j = 0; j < deploys[i].vflavors.size(); j++)
		{
			if (deploys[i].nums[j]) //not 0
			{
				result_file.append(deploys[i].vflavors[j]).append(" ");
				sprintf(s, "%d", deploys[i].nums[j]);
				result_file.append(s).append(" ");
			}
		}
		result_file.append("\n");
	}

	cout << result_file;
	write_result(result_file.c_str(), filename);


}

