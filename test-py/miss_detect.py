import sys
import os
import string


def check_dir(dir_name,out_file,exe_num):
	list_file = os.listdir(dir_name);
	for line in list_file:
		path_tmp = "%s\\%s"%(dir_name,line);
		if os.path.isdir(path_tmp):
			check_dir(path_tmp,out_file,exe_num);
		else:
			if path_tmp[-4:] == ".wav":
				print(path_tmp);
				cmd="%d.exe %s %s "%(exe_num,path_tmp, out_file);
				print(cmd)
				os.system(cmd);
	return 0;

if len(sys.argv) < 4:
    print("usage: %s num wav_dir out.log \n"%(sys.argv[0]));
    sys.exit(0)


exe_num=int(sys.argv[1]);
dir_root=sys.argv[2];
out_log=sys.argv[3]
fp_log=open(out_log,"w");
fp_log.close();
check_dir(dir_root,out_log,exe_num);

fp_log=open(out_log,"r");
c_all = 0;
c_e = 0;
lines=fp_log.readlines();
for line in lines:
	c_all += 1;
	vec=line.split("\t");
	if vec[1] == "1":
		c_e += 1;

fp_log.close()


fp_log=open(out_log,"a");
fp_log.write("all_num=%d\terr_num=%d\tratio=%.4f\n"%(c_all,c_e,float(c_e)/float(c_all)*100));
fp_log.close()






