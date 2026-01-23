from tro_org.benchmark.usage_profiler import profile_resources
import tro_org.benchmark.usage_profiler as profiler
import time

def counter(number):
    res = []
    for i in range(number):
        calc = (number*i)*(number *i)*(number *i)
        print(calc)
        res.append(calc)


def counter2(number):
    res = []
    for i in range(number):
        calc = (number*i)*(number *i)
        print(calc)
        res.append(calc)




if __name__ == '__main__':
    strt = time.time()
    rundef = profile_resources("test", "ztest_folder/test_usage.tsv")(counter)
    print(rundef(10000000))
    rundef2 = profile_resources("felix", "ztest_folder/felix_usage.tsv")(counter2)
    print(rundef2(10000000))
    now = time.time()
    print(now - strt)
