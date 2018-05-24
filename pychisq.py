#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from scipy import stats
from scipy.optimize import newton
from scipy.special import chdtri


def main():
    """
    Main Processing.
    """
    # intial setting --- variable definition.
    nc0 = 0.00001 # Do not change this variable.
    itr = 100000 # Do not change this variable.
    slesh_hold = 0.001 # Do not change this variable.

    w = 0.3 # Cohen's w : Can be changed.
    N = None # Sample size : Can be changed.
    df = 1 # Degrees of Freedom : Can be changed.
    sig_level = 0.05 # Significant Level: Can be changed.
    power = 0.9 # Power: Can be changed.

    var_list = [w, N, df, sig_level, power] # Put all variables in a list

    if(None_var_num_check(var_list = var_list) > 1):
        """
        Check whether or not an appropriate variables are introduced.
        """
        print('エラー：w, N, df, sig_level, power の中に' +
                      '数値の定まっていないものが２つ以上あります。' +
                      '１つに減らしてください')
        exit() # End Processing due to a completion of the porocess
    elif(None_var_num_check(var_list = var_list) == 0):
        print('エラー：求めるべき変数がないため、プログラムを終了します。')
        exit() # End Processing due to a completion of the porocess

    # check variables' varidity
    var_check(w = w, N = N, df = df, sig_level = sig_level, power = power)


    if(w is None):
        """
        Compute w
        """

        for i in range(itr):
            w0 = i/itr # Set initial w
            p = obj_function(w = w0, N = N, df = df
                             , sig_level = sig_level
                             , power = power, nc0 = nc0)
            if(abs(p) < slesh_hold):
                print('wの値は{0}です'.format(w0))
                exit() # End Processing due to a completion of the porocess

        if(p == False):
            for i in range(itr):
                w0 = - i/itr # Set initial w
                p = obj_function(w = w0, N = N, df = df
                                 , sig_level = sig_level
                                 , power = power, nc0 = nc0)
                if(abs(p) < slesh_hold):
                    print('wの値は{0}です'.format(w0))
                    exit() # End Processing due to a completion of the porocess

        if(p is None):
            print('申し訳ありません、wの値は求められませんでした。')
            exit() # End Processing due to a completion of the porocess


    if(N is None):
        """
        Compute N
        """
        for i in range(itr):
            N0 = i + 1 # Set initial N
            p = obj_function(w = w, N = N0, df = df
                             , sig_level = sig_level
                             , power = power, nc0 = nc0)
            if(abs(p) < slesh_hold * 2): # due to sensitivity multiplying 2.
                print('Nの値は{0}です'.format(N0))
                exit() # End Processing due to a completion of the porocess

        if(p == False):
            for i in range(itr):
                N0 = - i - 1
                p = obj_function(w = w, N = N0, df = df
                                 , sig_level = sig_level
                                 , power = power, nc0 = nc0)
                if(abs(p) < slesh_hold * 2): # due to sensitivity multiplying 2.
                    print('Nの値は{0}です'.format(N0))
                    exit() # End Processing due to a completion of the porocess
    
        if(i == itr - 1):
            print('申し訳ありません、Nの値は求められませんでした。')
            exit() # End Processing due to an Error

# The followings are functions.

def var_check(w = None, N = None, df = None, sig_level = None, power = None):
    """
    Check wheter or not appropriate numbers are introduced.
    """   
    if(w is not None and isinstance(w, float)):
        if(0 >= w or w >= 1):
            print('エラー：wは0～1の間の実数である必要があります。' +
                          'プログラムを終了します')
            exit() # End Processing due to an Error
    elif(w is None):
        print('wは入力されていません。wを求めます')

    if(N is not None and N <= 0):
        print('エラー：Nは1以上の整数である必要があります。' +
                           'プログラムを修理します')
        exit() # End Processing due to an Error
    elif(N is None):
        print('Nは入力されていません。Nを求めます')

    if(df is not None and df <= 0):
        print('dfには1以上の数字を入れて下さい')
    elif(df is None):
        print('エラー：dfは入力されていません。プログラムを終了します')
        exit() # End Processing due to an Error

    if(sig_level is not None and isinstance(sig_level, float)):
        if(0 > sig_level or sig_level > 1):
            print('エラー：sig_levelは0～1の間のである必要があります' +
                       'プログラムを終了します。')
            exit() # End Processing due to an Error
    elif(sig_level is None):
        print('エラー：sig_levelは入力されていません。プログラムを終了します')
        exit() # End Processing due to an Error

    if(power is not None and isinstance(power, float)):
        if(0 > power or power > 1):
            print('エラー：powerは0～1の間の実数である必要があります' +
                       'プログラムを終了します。')
            exit() # End Processing due to an Error
    else:
        print('powerは入力されていません。powerを求めます')


def None_var_num_check(var_list = None):
    """
    Check whether or not there are missing values more than two in var_list.
    """   
    a = 0
    for i in var_list:
        if i is None:
            a += 1
    return a 


def py_pchisq(q = None, df = None, ncp = None
              , lower_tail = True, log_p = False, nc0 = None):

    if(ncp > 0):
        s = 1 - stats.ncx2.cdf(q, df, nc = ncp)
    else:
        s = 1 - stats.ncx2.cdf(q, df, nc = nc0)
    return s


def power_compute(w = None, N = None, df = None
                  , sig_level = None, power = None, nc0 = None):
    """
    Compute power with other variables.
    """   
    k = None # set k
    if(power is None):
        k = chdtri(df, sig_level) # Compute percentile
    return py_pchisq(q = k, df = df, ncp = N*w**2, nc0 = nc0)

def obj_function(w = None, N = None, df = None
                 , sig_level = None, power = None, nc0 = None):
    """
    Objective function which is used for the w and N computation.
    """   
    k = chdtri(df, sig_level) # Set cumulative distribution function.
    k = py_pchisq(q = k, df = df, ncp = N*w**2, nc0 = nc0) - power
    return k




if __name__ == '__main__':
    main()
