#! /usr/bin/python3

import os
import os.path

class DensityParser:
    density_names = ['sparse', 'intermediate', 'dense']
    @classmethod
    def __classify_density_lh(cls, densestr):
        density = 1
        if densestr[1:].startswith('-1'):           # low value is 1
            density = 0
        elif densestr.startswith('1000-'):   # high value is 1000
            density = 2
        return density
    @classmethod
    def __classify_density(cls, densestr):
        for i, dnstr in enumerate(cls.density_names):
            if densestr == dnstr:
                return i
        return 1

    @classmethod
    def parse_command_line(cls, time_cmd):
        ''' parse density result. e
        time_cmd = 
        Command being timed: "gt_mpi_gather -j 1-qc_pos_1000_1-1.json --print-calls -s 5000"
        > e-9: 16-qc_dense_1000_678-1000.json
        return [2, 10, 3000] => dense, 10 pos, seg-size 3000
        '''
        cmds = time_cmd.split()
        token = None
        if len(cmds) >= 6:
            for i, cmd in enumerate(cmds):
                if cmd.endswith('gt_mpi_gather'):
        # if len(cmds) == 6 and cmds[0].endswith('gt_mpi_gather'):  # only query
            # print("INFO: +++cmd: %s" % time_cmd)
                    qc_name = os.path.basename(cmds[i+2]).split('_')
                    density = cls.__classify_density_lh(qc_name[-1]) if i == 0 and len(cmds) < 8 else cls.__classify_density(qc_name[1])
                    #[density, num_pos, seg_size]
                    token = tuple((density, int(qc_name[-2]), int(cmds[i+5]), cmds[i+3]))  
                    #print("INFO: ---parse_density return: %s" % cls.to_report_str(token))
                    break
        # else:
        #     print('INFO skipped non-"gt_mpi_gather", cmd %s' % time_cmd)
        return token

    @classmethod
    def to_report_str(cls, t):
        ''' from [2, 10, 3000] => dense 10 pos with seg_size 3000 '''
        return "%s %d pos with seg_size %d" % (cls.density_names[t[0]], t[1], t[2])
