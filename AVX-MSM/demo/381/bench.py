import pandas as pd
import subprocess

def deal_result(bench_res: str):
    t = bench_res.strip()
    t = t.split("\n")
    res = []
    
    for item in t[:2]:
        res.append(float(item.split("=")[1]))
    res.append(res[-2]/res[-1])
    if t[2] != "YES":
        print("[FAIL] Result Error")
        exit(1)
    return res

def bench_pip_ifma(num=10, wsize=4, nbench=1):
    try:
        subprocess.check_call(
            'make clean',
            shell=True,
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
        )
        subprocess.run(
            f'make ifma CDEFINES="-DWNUM={num} -DWMBITS={wsize} -DNBENCHS={nbench}"', 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True,
        )
    except:
        print("[FAIL] Compile Error")
        return
    
    try:
        
        result = subprocess.run(
            ['build/test_pip_ifma'], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True,
        )
        # print(result.stdout)
        
        tmp = deal_result(result.stdout)
        tmp.insert(0, wsize)
        tmp.insert(0, num)
        # print(tmp)

        return tmp
    except KeyboardInterrupt:
        print("[Stop] KeyboardInterrupt")
        exit(1)
    except:
        print("[FAIL] Run Error")
        exit(1)

def bench_pip_thread(num=10, wsize=4, nbench=1):
    try:
        subprocess.check_call(
            'make clean',
            shell=True,
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
        )
        subprocess.run(
            f'make thread CDEFINES="-DWNUM={num} -DWMBITS={wsize} -DNBENCHS={nbench}"', 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True,
        )
    except:
        print("[FAIL] Compile Error")
        return
    
    try:
        
        result = subprocess.run(
            ['build/test_pip_threads'], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True,
        )
        # print(result.stdout)
        
        tmp = deal_result(result.stdout)
        tmp.insert(0, wsize)
        tmp.insert(0, num)
        # print(tmp)

        return tmp
    except KeyboardInterrupt:
        print("[Stop] KeyboardInterrupt")
        exit(1)
    except:
        print("[FAIL] Run Error")
        exit(1)

def bench_pair_ifma(num=10, wsize=4, nbench=1):
    try:
        subprocess.check_call(
            'make clean',
            shell=True,
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
        )
        subprocess.run(
            f'make pair_ifma CDEFINES="-DWNUM={num} -DWMBITS={wsize} -DNBENCHS={nbench}"', 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True,
        )
    except:
        print("[FAIL] Compile Error")
        return
    
    try:
        
        result = subprocess.run(
            ['build/test_pair_ifma'], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True,
        )
        # print(result.stdout)
        
        tmp = deal_result(result.stdout)
        tmp.insert(0, wsize)
        tmp.insert(0, num)
        # print(tmp)

        return tmp
    except KeyboardInterrupt:
        print("[Stop] KeyboardInterrupt")
        exit(1)
    except:
        print("[FAIL] Run Error")
        exit(1)

def bench_pair_thread(num=10, wsize=4, nbench=1):
    try:
        subprocess.check_call(
            'make clean',
            shell=True,
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
        )
        subprocess.run(
            f'make pair_thread CDEFINES="-DWNUM={num} -DWMBITS={wsize} -DNBENCHS={nbench}"', 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True,
        )
    except:
        print("[FAIL] Compile Error")
        return
    
    try:
        
        result = subprocess.run(
            ['build/test_pair_threads'], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True,
        )
        # print(result.stdout)
        
        tmp = deal_result(result.stdout)
        tmp.insert(0, wsize)
        tmp.insert(0, num)
        # print(tmp)

        return tmp
    except KeyboardInterrupt:
        print("[Stop] KeyboardInterrupt")
        exit(1)
    except:
        print("[FAIL] Run Error")
        exit(1)



if __name__ == '__main__':
    df = pd.DataFrame(columns=["NUM", "WSIZE", "OLD", "NEW", "SPEEDUP"])
    df_best = pd.DataFrame(columns=["NUM", "WSIZE", "OLD", "NEW", "SPEEDUP"])

    num_max = 25
    wsize_gap1 = 1
    wsize_gap2 = 9
    nbench = 1

    for num in range(15, 16):
        max_time = 10000000000000000
        best = []
        for wsize in range(num - wsize_gap2, num - wsize_gap1):
            # print(f"NUM={num}, WSIZE={wsize}")
            run_time = bench_pip_thread(num, wsize, nbench)
            print(run_time)
            df.loc[len(df)] = run_time

            if run_time[3] < max_time:
                max_time = run_time[3]
                best = run_time
        df_best.loc[len(df_best)] = best
        print(f"Best: {best}")
    
    df_best.to_csv("best_thread.csv")
    df.to_csv("result_thread.csv")