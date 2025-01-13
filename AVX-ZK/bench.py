import pandas as pd
import subprocess

def deal_result(bench_res):
    if 'The verification result is: PASS' not in bench_res:
        return False, None
    t = bench_res.strip()
    t = t.split("\n")
    
    tag = '(leave) Call to r1cs_gg_ppzksnark_prover'
    for s in t:
        if tag in s:
            return True, s[s.index('['):]
    return False, None

def bench_zksnark(size1=1000, size2=100):
    
    try:
        
        result = subprocess.run(
            [f'./build/libsnark/profile_r1cs_gg_ppzksnark {size1} {size2} bytes'], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True,
            text=True
        )
        # print(result.stdout)
        
        flag, res =  deal_result(result.stdout)
        tmp = [flag, res]

        return tmp
    except KeyboardInterrupt:
        print("[Stop] KeyboardInterrupt")
        exit(1)



if __name__ == '__main__':
    df = pd.DataFrame(columns=["SIZE", "Correct", "Time"])
    # df_best = pd.DataFrame(columns=["NUM", "WSIZE", "OLD", "NEW", "SPEEDUP"])for
    # bench_zksnark(1000, 100)
    for num in range(23, 25):
        input1 = (1 << num)
        input2 = (1 << (num-2))
        run_time = bench_zksnark(input1, input2)
        run_time.insert(0, num)
        print(run_time)
        df.loc[len(df)] = run_time

    df.to_csv("result.csv")