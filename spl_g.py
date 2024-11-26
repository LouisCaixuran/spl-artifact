import re
import csv


PATH="./"

# Function to process the input file
def process_file(filename):
    with open(filename, 'r') as file:
        raw_lines = file.readlines()

    lines=[]
    for i in range(len(raw_lines)):
        if raw_lines[i].startswith('linear_scan_r'):
            lines.append(raw_lines[i-1])
            lines.append(raw_lines[i])

    data = []
    for i in range(0, len(lines), 2):

        cfg_line = lines[i].strip()
        time_line = lines[i+1].strip()

        # Extracting register values
        cfg_size=-1
        try:
            cfg_size=int(re.search(r'cfg size = (\d+)', cfg_line).group(1))
        except:
            print(cfg_line)

        registers = {
            'linear_scan_r': int(re.search(r'linear_scan_r: (-?\d+)', time_line).group(1)),
            'graph_color_r': int(re.search(r'graph_color_r: (-?\d+)', time_line).group(1)),
            'spl_r': int(re.search(r'spl_r: (-?\d+)', time_line).group(1)),
            'treedec_r': int(re.search(r'treedec_r: (-?\d+)', time_line).group(1)),
            'pathdec_r': int(re.search(r'pathdec_r: (-?\d+)', time_line).group(1)),
            'puzzle_r': int(re.search(r'puzzle_r: (-?\d+)', time_line).group(1)),
        }

        # Extracting time values
        times = {
            'linear_scan_time': int(re.search(r'linear_scan_time: (\d+)', time_line).group(1)),
            'graph_color_time': int(re.search(r'graph_color_time: (\d+)', time_line).group(1)),
            'spl_time': int(re.search(r'spl_time: (\d+)', time_line).group(1)),
            'treedec_time': int(re.search(r'treedec_time: (\d+)', time_line).group(1)),
            'pathdec_time': int(re.search(r'pathdec_time: (\d+)', time_line).group(1)),
            'puzzle_time': int(re.search(r'puzzle_time: (\d+)', time_line).group(1)),
        }

        # Append the data as a dictionary
        data.append({"cfg_size": cfg_size , **registers, **times})

    return data


def print_info(data):
    print("total case:" + str(len(data)))
    print("invalid spl_case:" + str(len([d for d in data if d['spl_r'] == -1])))
    print("invalid linear_scan_case:" + str(len([d for d in data if d['linear_scan_r'] == -1])))
    print("invalid graph_color_case:" + str(len([d for d in data if d['graph_color_r'] == -1])))
    print("invalid spl_case:" + str(len([d for d in data if d['spl_r'] == -1])))
    print("invalid treedec_case:" + str(len([d for d in data if d['treedec_r'] == -1])))
    print("invalid pathdec_case:" + str(len([d for d in data if d['pathdec_r'] == -1])))
    print("invalid puzzle_case:" + str(len([d for d in data if d['puzzle_r'] == -1])))

    #write the original data into csv
    with open(PATH+'data.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['cfg_size','opt_r', 'linear_scan_r', 'linear_scan_time', 'graph_color_r', 'graph_color_time', 'spl_r', 'spl_time', 'treedec_r', 'treedec_time', 'pathdec_r', 'pathdec_time', 'puzzle_r', 'puzzle_time'])
        for d in data:
            if d['graph_color_r'] == -1 :
                d['graph_color_time']=60000000
            if d['spl_r'] == -1 :
                d['spl_time']=60000000
            if d['treedec_r'] == -1 :
                d['treedec_time']=60000000
            if d['pathdec_r'] == -1 :
                d['pathdec_time']=60000000
            opt_r=max(d['treedec_r'],d['spl_r'],d['pathdec_r'],d['graph_color_r'])
            if opt_r==-1 or opt_r==0:
                continue
            writer.writerow([d['cfg_size'],opt_r, d['linear_scan_r'], d['linear_scan_time'], d['graph_color_r'], d['graph_color_time'], d['spl_r'], d['spl_time'], d['treedec_r'], d['treedec_time'], d['pathdec_r'], d['pathdec_time'], d['puzzle_r'], d['puzzle_time']])


    #write spl vs graphC
    with open (PATH+'spl_vs_graphC.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['opt_r', 'spl_time', 'graph_color_time','graph/spl','spl_invalid','graph_color_invalid'])
        for d in data:
            opt_r=max(d['spl_r'],d['graph_color_r'])
            if opt_r==-1 or opt_r==0:
                continue
            division=None
            spl_invalid=None
            graph_color_invalid=None
            if d['spl_r'] == -1:
                spl_invalid=0.1
                division=0.1
            elif d['graph_color_r'] == -1:
                graph_color_invalid=100000
                division=100000
            else:
                division=d['graph_color_time']/d['spl_time']
            writer.writerow([opt_r, d['spl_time'], d['graph_color_time'],division,spl_invalid,graph_color_invalid])

    #write treedec vs spl
    with open (PATH+'spl_vs_tree.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['opt_r', 'spl_time', 'graph_color_time','tree/spl','spl_invalid','tree_invalid'])
        for d in data:
            opt_r=max(d['spl_r'],d['treedec_r'])
            if opt_r==-1 or opt_r==0:
                continue
            division=None
            spl_invalid=None
            graph_color_invalid=None
            if d['spl_r'] == -1:
                spl_invalid=0.1
                division=0.1
            elif d['treedec_r'] == -1:
                graph_color_invalid=1500000
                division=1500000
            else:
                division=d['treedec_time']/d['spl_time']
            writer.writerow([opt_r, d['spl_time'], d['treedec_time'],division,spl_invalid,graph_color_invalid])


    #write pathdec vs spl
    with open (PATH+'spl_vs_path.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['opt_r', 'spl_time', 'graph_color_time','path/spl','spl_invalid','path_invalid'])
        for d in data:
            opt_r=max(d['spl_r'],d['pathdec_r'])
            if opt_r==-1 or opt_r==0:
                continue
            division=None
            spl_invalid=None
            graph_color_invalid=None
            if d['spl_r'] == -1:
                spl_invalid=0.1
                division=0.1
            elif d['pathdec_r'] == -1:
                graph_color_invalid=150000
                division=150000
            else:
                division=d['pathdec_time']/d['spl_time']
            writer.writerow([opt_r, d['spl_time'], d['pathdec_time'],division,spl_invalid,graph_color_invalid])


    #collect the optimal register value
    r=[0 for i in range(30)]
    for i in data:
        opt=max([i['treedec_r'],i['spl_r'],i['pathdec_r'],i['graph_color_r']])
        if opt!=-1:
            r[opt]+=1
    data_without_0=[d for d in data if max(d['treedec_r'],d['spl_r'],d['pathdec_r'],d['graph_color_r'])!=0]
    print("all case without 0: " + str(len(data_without_0)))
    with open(PATH+'optimal_register.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['register', 'count'])
        for i in range(30):
            writer.writerow([i, r[i]])


    data_all_valid=[d for d in data if d['linear_scan_r'] != -1 and d['graph_color_r'] != -1 and d['spl_r'] != -1 and d['treedec_r'] != -1 and d['pathdec_r'] != -1 and d['puzzle_r'] != -1]
    print("all valid case:" + str(len(data_all_valid)))
    print("average linear_scan_time:" + str(sum([d['linear_scan_time'] for d in data_all_valid])/len(data_all_valid)))
    print("average graph_color_time:" + str(sum([d['graph_color_time'] for d in data_all_valid])/len(data_all_valid)))
    print("average spl_time:" + str(sum([d['spl_time'] for d in data_all_valid])/len(data_all_valid)))
    print("average treedec_time:" + str(sum([d['treedec_time'] for d in data_all_valid])/len(data_all_valid)))
    print("average pathdec_time:" + str(sum([d['pathdec_time'] for d in data_all_valid])/len(data_all_valid)))
    print("average puzzle_time:" + str(sum([d['puzzle_time'] for d in data_all_valid])/len(data_all_valid)))





    data_fail_linear_scan=[d for d in data if (d['linear_scan_r']>d['treedec_r'] and d['treedec_r']!=-1) or (d['linear_scan_r']>d['spl_r'] and d['spl_r']!=-1)]
    data_fail_puzzle=[ d for d in data if (d['puzzle_r']>d['pathdec_r'] and d['pathdec_r']!=-1) or (d['puzzle_r']>d['spl_r'] and d['spl_r']!=-1)]
    print("fail linear_scan register:" + str(len(data_fail_linear_scan)))
    print("fail puzzle register:" + str(len(data_fail_puzzle)))
    with open(PATH+'fail_linear_scan.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['cfg', 'linear_scan_r','optimal_r','differece'])
        for d in data_fail_linear_scan:
            writer.writerow([d['cfg_size'], d['linear_scan_r'], max(d['treedec_r'],d['spl_r']), d['linear_scan_r']-max(d['treedec_r'],d['spl_r'])])


    data_all_valid_without_0=[d for d in data_all_valid if d['treedec_r']!=0]

    print("all valid case without 0"+ str(len(data_all_valid_without_0)))
    #print avertege time without 0
    print("average linear_scan_time without 0:" + str(sum([d['linear_scan_time'] for d in data_all_valid_without_0])/len(data_all_valid_without_0)))
    print("average graph_color_time without 0:" + str(sum([d['graph_color_time'] for d in data_all_valid_without_0])/len(data_all_valid_without_0)))
    print("average spl_time without 0:" + str(sum([d['spl_time'] for d in data_all_valid_without_0])/len(data_all_valid_without_0)))
    print("average treedec_time without 0:" + str(sum([d['treedec_time'] for d in data_all_valid_without_0])/len(data_all_valid_without_0)))
    print("average pathdec_time without 0:" + str(sum([d['pathdec_time'] for d in data_all_valid_without_0])/len(data_all_valid_without_0)))
    print("average puzzle_time without 0:" + str(sum([d['puzzle_time'] for d in data_all_valid_without_0])/len(data_all_valid_without_0)))

    #write into file
    with open(PATH+'times_without_0.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['cfg_size', 'opt_r','graph_color_time','spl_time','treedec_time','pathdec_time'])
        for d in data_all_valid_without_0:
            writer.writerow([d['cfg_size'],d['treedec_r'], d['graph_color_time'],d['spl_time'],d['treedec_time'],d['pathdec_time']])

    # count the number of registers:
    opt_r=[0 for i in range(21)]
    puzzle_r=[0 for i in range(21)]
    linear_r=[0 for i in range(21)]
    diff_register_linear=[[0 for i in range(21)] for i in range(21)]
    diff_register_puzzle=[[0 for i in range(21)] for i in range(21)]


    for d in data:
        opt=max(d['treedec_r'],d['spl_r'],d['pathdec_r'],d['graph_color_r'])
        if opt==-1 or d['linear_scan_r']==0:
            continue
        opt_r[opt-1]+=1
        if d['puzzle_r']<=20 and d['puzzle_r']!=-1:
            puzzle_r[d['puzzle_r']-1]+=1
            diff_register_puzzle[opt-1][d['puzzle_r']-1]+=1
        else:
            puzzle_r[20]+=1
            diff_register_puzzle[opt-1][20]+=1
        if d['linear_scan_r']<=20:
            linear_r[d['linear_scan_r']-1]+=1
            diff_register_linear[opt-1][d['linear_scan_r']-1]+=1
        else:
            linear_r[20]+=1
            diff_register_linear[opt-1][20]+=1

    with open(PATH+'register_count_c.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['register', 'opt_r','puzzle_r','linear_r'])
        for i in range(21):
            writer.writerow([i+1, opt_r[i],puzzle_r[i],linear_r[i]])

    with open(PATH+'diff_register_linear.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow([(i+1) for i in range(20)])
        for i in range(20):
            writer.writerow(diff_register_linear[i])

    with open(PATH+'diff_register_puzzle.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow([(i+1) for i in range(20)])
        for i in range(20):
            writer.writerow(diff_register_puzzle[i])

    with open(PATH+'fail_puzzle.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['cfg', 'puzzle_r','optimal_r','differece'])
        for d in data_fail_puzzle:
            opt_r=max(d['pathdec_r'],d['spl_r'])
            writer.writerow([d['cfg_size'], d['puzzle_r'], opt_r, d['puzzle_r']-opt_r])
    return data_all_valid

def collect_diff(d1,d2):
    with open(d1, 'r') as file:
        lines_linear = file.readlines()
    with open(d2, 'r') as file:
        lines_puzzle= file.readlines()

    d_linear=[0 for i in range(8)]
    d_puzzle=[0 for i in range(8)]

    #read diff from linear_scan
    for i in range(1,len(lines_linear)):
        diff=int(lines_linear[i].split(',')[3])
        d_linear[diff-1]+=1

    #read diff from puzzle
    for i in range(1,len(lines_puzzle)):
        diff=int(lines_puzzle[i].split(',')[3])
        d_puzzle[diff-1]+=1


    #write into csv file
    with open(PATH+'diff_register.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['diff', 'count_linear', 'count_puzzle'])
        for i in range(8):
            writer.writerow([i+1, d_linear[i],d_puzzle[i]])





# Main execution
if __name__ == "__main__":
    filename = PATH+"output.txt"  # Replace with your actual file name
    data = process_file(filename)
    data_all_valid=print_info(data)

    collect_diff(PATH+'fail_linear_scan.csv',PATH+'fail_puzzle.csv')
    #write times and cfg size to csv file
    with open(PATH+'times.csv', mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['cfg_size', 'graph_color_time','spl_time','treedec_time','pathdec_time'])
        for d in data_all_valid:
            writer.writerow([d['cfg_size'], d['graph_color_time'],d['spl_time'],d['treedec_time'],d['pathdec_time']])


