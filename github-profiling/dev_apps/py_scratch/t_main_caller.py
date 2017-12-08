#! /usr/bin/python
from t_caller import f, TestCaller
if __name__ == "__main__":
    f()
    print("---TestCaller")
    TestCaller().f()
