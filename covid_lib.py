# -*- coding: utf-8 -*-
"""
日付:
    2020年07月26日
    2021年01月21日
    2021年02月06日
		2021年08月15日
ファイル名:
    covid_lib.py
内容:
    新型コロナウイルスの新規患者数と(簡易・日毎)実効再生産数を計算するプログラム
入力:
    北海道、東京、横浜、大阪、福岡のCSVファイル
    NHKのCSVファイル
出力:
    graph.png
"""
import datetime
import datetime as dt
import decimal
import sys

import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
import scipy.optimize
from matplotlib.figure import Figure
from scipy import optimize

"""最小二乗法で移動平均をとるための関数群"""


# (15-1)二次関数、単一変数版
def parabola(z, a, b, c):
    return a * z ** 2 + b * z + c


# (15-2)直線、単一変数版
def linear(z, a, b):
    return a * z + b


# (15-3)10を底とする対数関数、負の場合はNaNを返す、単一変数版
def dB10_float(z):
    return np.log10(z) if z > 0 else np.NaN


# (15-4)10を底とする冪乗関数、単一変数版
def power10_float(z):
    if z == z:
        # print (z)
        return 10 ** (z)
    else:
        return np.NaN
    # return 10 ** (z / 10)


dB10_df = np.frompyfunc(dB10_float, 1, 1)  # (15-5)10を底とする対数関数、負の場合はNaNを返す、numpy版

power10_df = np.frompyfunc(power10_float, 1, 1)  # (15-6)10を底とする冪乗関数、numpy版


# (15-7)感染者数の二次関数カーブフィットするオブジェクト1
def fit_by_parabola(x_se, y_se, extend_length):
    param_se, pcov_se = scipy.optimize.curve_fit(parabola, x_se, y_se)  # カーブのパラメータをscipyの最適化関数で計算
    return fit_se  # , param_se, perr


# (15-8)感染者数のカーブを区間t0毎に分け、区間毎に最小自乗近似して、各区間推定カーブの右端の値のみを取って全区間の推定値を計算するオブジェクト
def fit_by_linear_all_range(x, y, t_segment, t_extension):
    wn = y.copy()  # expanded x
    wn = pd.Series(np.zeros(len(y)))
    zn = pd.Series(np.zeros(t_segment))
    x = x.interpolate()  # 補間
    y = y.interpolate()
    x = x.fillna(0)  # ０補間
    y = y.fillna(0)
    xn = pd.Series()
    yn = pd.Series()
    for n in range(t_segment - 1, len(x)):
        xn = x[n - t_segment + 1:n].copy()  # 区間t0コピー
        yn = y[n - t_segment + 1:n].copy()
        param_se, pcov_se = scipy.optimize.curve_fit(linear, xn, yn)  # カーブのパラメータをscipyの最適化関数で計算
        # perr = np.sqrt(np.diag(pcov_se))  # 誤差の計算
        x_extension_se = pd.Series(range(xn.iloc[0], xn.iloc[0] + len(xn) + t_extension + 1))
        zn = linear(x_extension_se, param_se[0], param_se[1])  # Curve に計算したパラメータを入れて近似曲線を計算する
        wn.iloc[n] = zn.iloc[-1]  # 最後のデータを取る
    return wn  # 推定した曲線


# (15-9)感染者数のカーブを区間t0毎に分け、区間毎に最小自乗近似して、各区間推定カーブの右端の値のみを取って全区間の推定値を計算するオブジェクト
def fit_by_parabola_all_range(x, y, t_segment, t_extension):
    wn = y.copy()  # expanded x
    wn = pd.Series(np.zeros(len(y)))
    zn = pd.Series(np.zeros(t_segment))
    x = x.interpolate()  # 補間
    y = y.interpolate()
    x = x.fillna(0)  # ０補間
    y = y.fillna(0)
    xn = pd.Series()
    yn = pd.Series()
    for n in range(t_segment - 1, len(x)):
        xn = x[n - t_segment + 1:n].copy()  # 区間t0コピー
        yn = y[n - t_segment + 1:n].copy()
        param_se, pcov_se = scipy.optimize.curve_fit(parabola, xn, yn)  # カーブのパラメータをscipyの最適化関数で計算
        # perr = np.sqrt(np.diag(pcov_se))  # 誤差の計算
        x_extension_se = pd.Series(range(xn.iloc[0], xn.iloc[0] + len(xn) + t_extension + 1))
        zn = parabola(x_extension_se, param_se[0], param_se[1], param_se[2])  # Curve に計算したパラメータを入れて近似曲線を計算する
        wn.iloc[n] = zn.iloc[-1]  # 最後のデータを取る
    return wn  # 推定した曲線


# (15-10)感染者数のカーブを区間t0毎に分け、区間毎に最小自乗近似して、各区間推定カーブの右端の値のみを取って全区間の推定値を計算するオブジェクト
def median_all_range(y, t0):
    wn = y.copy()  # expanded x
    wn = pd.Series(np.zeros(len(y)))
    for n in range(t0 - 1, len(y)):
        wn.iloc[n] = np.median(y[n - t0 + 1:n])  # 区間のメディアンを計算
    return wn  # メディアンした曲線


# (15-11)感染者数のカーブを区間t0毎に分け、区間毎に最小自乗近似して、各区間推定カーブの右端の値のみを取って全区間の推定値を計算するオブジェクト
def max_all_range(y, t0):
    wn = y.copy()  # expanded x
    wn = pd.Series(np.zeros(len(y)))
    for n in range(t0 - 1, len(y)):
        wn.iloc[n] = np.max(y[n - t0 + 1:n])  # 区間のメディアンを計算
    return wn  # メディアンした曲線


# (15-12)感染者数のカーブを区間t0毎に分け、区間毎に最小自乗近似して、各区間推定カーブの右端の値のみを取って全区間の推定値を計算するオブジェクト
def min_all_range(y, t0):
    wn = y.copy()  # expanded x
    wn = pd.Series(np.zeros(len(y)))
    for n in range(t0 - 1, len(y)):
        wn.iloc[n] = np.min(y[n - t0 + 1:n])  # 区間のメディアンを計算
    return wn  # メディアンした曲線


"""2つのベクトルのノルム２を最小化する関数"""


def find_x_min(data_df, model_df):
    # import scipy.optimize.brent
    def f(x):
        return np.sqrt((np.sum(data_df - (model_df + x)) ** 2))

    x_min = optimize.brent(f)  # It actually converges in 9 iterations!
    return x_min


"""メイン関数"""


class Covid19(object):
    def __init__(self):  # オブジェクト
        # 初期値
        self.df_df = pd.DataFrame()
        self.encoding = 'utf-8'
        self.area = '東京'
        self.title = 'Tokyo'
        # self.date_key = '公表_年月日'
        # self.date_key = '月日'
        self.date_key = '各地の感染者数_1日ごとの発表数'
        self.file_name = 'nhk_news_covid19_prefectures_daily_data.csv'
        self.New_cases = '各地の感染者数_1日ごとの発表数'
        self.Total_cases = '各地の感染者数_累計'
        self.pref = '県名'
        self.dir = './NHK全国データ/'
        self.date = '2021/03/20'

        # データ構造
        self.col = ['No', 'date', 'New_cases', 'Rt', 'i14', 'New_cases_rolling', 'Total_cases']
        self.cvd_df = pd.DataFrame(columns=self.col)

        # self.make_text_for_twitter
        self.prefectures_se = [
            '北海道', '青森県', '岩手県', '宮城県', '秋田県', '山形県', '福島県',
            '茨城県', '栃木県', '群馬県', '埼玉県', '千葉県', '東京都', '神奈川県',
            '新潟県', '富山県', '石川県', '福井県', '山梨県', '長野県',
            '岐阜県', '静岡県', '愛知県', '三重県', '京都府',
            '滋賀県', '大阪府', '兵庫県', '奈良県', '和歌山県',
            '鳥取県', '島根県', '岡山県', '広島県', '山口県',
            '徳島県', '香川県', '愛媛県', '高知県',
            '福岡県', '佐賀県', '長崎県', '熊本県', '大分県', '宮崎県', '鹿児島県', '沖縄県']

        self.emp_se = ['大阪府', '東京都', '兵庫県',

                       '埼玉県', '沖縄県', '宮城県', '神奈川県',  # 5-8
                       '愛知県',

                       '京都府', '奈良県', '千葉県', '北海道',  # 9-12
                       '長野県', '茨城県', '福岡県', '和歌山県',  # 13-16

                       # '福岡県',  '群馬県', '福島県',
                       # '栃木県', '滋賀県', '新潟県', '静岡県', '愛媛県', '山形県',
                       '広島県'
                       ]  # 新規患者数順
        # print(self.emp_se)
        # for excel output
        self.prefecture_names_se = []
        self.items_se = ['実効再生産数', 'ソートずみ実効再生産数県別順位', 'ソートずみ実効再生産数', '新規感染者数', 'ソートずみ新規感染者数県別順位']

        self.new_cases_all_df = pd.DataFrame(dtype=object)
        self.rt_all_df = pd.DataFrame(dtype=object)
        self.url_IRYO = "https://www.kantei.go.jp/jp/content/IRYO-vaccination_data.xlsx"
        self.url_KOREI = "https://www.kantei.go.jp/jp/content/KOREI-vaccination_data.xlsx"
        self.jinkou = 127138033  # CIO【総計】令和2年１月１日住民基本台帳年齢階級別人口（市区町村別）
        self.jinkou_64 = 91651156  # 64歳以下
        self.jinkou_65 = 35486813  # 65歳以上
        self.vc_df = pd.DataFrame()
        print('Covid19 class constructor')

    # (0)感染者数のcsv生データを読むオブジェクト
    def read_csv_file(self):
        # 2020-04-01形式、エクセル禁止2020/4/1はダメ
        self.df_df = pd.DataFrame(dtype=object)
        try:
            self.df_df = pd.read_csv(self.file_name, encoding=self.encoding)  # 2020-04-01形式、エクセル禁止、2020/4/1はダメ
            # print(self.df_df.head())
        except IOError:
            print('ごめんなさい、', self.file_name, 'が読めないので終了しますね')
            exit('異常終了')
        return self

    # (1)感染者数のexcel生データを読む
    def read_excel_file(self):  # 大阪用オブジェクト
        self.df_df = pd.DataFrame(dtype=object)
        try:
            self.df_df = pd.read_excel(self.file_name, header=1)  # 大阪府はエクセル、最初の1行を読み飛ばす
        except IOError:
            print('ごめんなさい、', self.file_name, 'が読めないので終了しますね')
            exit('異常終了')
            return self

    # (2) # １日の感染数から、累積感染者数を出す
    def count(self):  # １日の感染数から、累積感染者数を出す#県別データオブジェクト
        self.cvd_df['New_cases'] = self.df_df[self.date_key].value_counts()  # 日付ごとの新規感染者数数をカウント
        self.cvd_df = self.cvd_df.sort_index(ascending=True)  # 日付でソート
        self.cvd_df['date'] = self.cvd_df.index  # 日付をいれる
        # print('count')
        # print(self.cvd_df)
        return self

    # (3) # 日付ごとの感染者数をcvd_dfへコピー
    def copy(self):  # 日付ごとの感染者数をcvd_dfへコピー＃NHKデータオブジェクト
        # print('copy')
        self.cvd_df['New_cases'] = self.df_df[self.date_key].value_counts()  # 日付ごとの新規感染者数数をカウント
        self.cvd_df['date'] = self.df_df[self.date_key].replace('/', '-', regex=True)  # 2020/01/01 から2020-01-01　へ変換
        self.cvd_df = self.cvd_df.sort_index(ascending=True)  # 日付でソート
        self.cvd_df.set_index('date')  # 日付をいれる
        # print('copy')
        # print(self.cvd_df)
        return self

    # (4)感染者数の生データをCUMSUM,Rt処理オブジェクト
    def compute_rt(self):  # 累積感染者数から７日平均,Rtを出す
        # print('compute_rt\n')
        self.cvd_df['No'] = range(len(self.cvd_df['New_cases']))  # 連番を作る
        self.cvd_df['Total_cases'] = self.cvd_df['New_cases'].cumsum()  # 累積感染者数を計算
        self.cvd_df['New_cases_rolling'] = self.cvd_df['New_cases'].rolling(7, min_periods=7).sum() / 7.0  # ７日平均をとる

        self.sw = 1
        if self.sw == 0:
            self.cvd_df['New_cases_rolling'] = self.cvd_df['New_cases'].rolling(7, min_periods=7).sum() / 7.0  # ７日平均をとる
            self.cvd_df['Rt'] = (self.cvd_df['New_cases_rolling'].pct_change(periods=7) + 1) ** (5 / 7)
        elif self.sw == 1:  # dB領域で、Rtと平均的なトレースを計算する
            self.cvd_df['New_cases_dB'] = dB10_df(self.cvd_df['New_cases'])  # dBへ変換
            self.cvd_df['Rt_1_dB'] = (self.cvd_df['New_cases_dB'].diff(7)) * (1 / 7)  # dB領域で１日勾配を計算
            self.cvd_df['Rt_5_dB'] = (self.cvd_df['New_cases_dB'].diff(7)) * (5 / 7)  # dB領域でRtを計算
            self.cvd_df['New_cases_rolling_dB'] = self.cvd_df['Rt_1_dB'].cumsum()  # １日勾配を積分して平均を計算
            x_min = find_x_min(self.cvd_df['New_cases_dB'], self.cvd_df['New_cases_rolling_dB'])  # DC補正値を計算、Norm2
            self.cvd_df['New_cases_rolling_dB'] += x_min  # DC補正値を加算
            average_days = 1  # アベレージの日数をセット
            self.cvd_df['Rt_5_dB'] = self.cvd_df['Rt_5_dB'].rolling(average_days,
                                                                    min_periods=average_days).sum() / average_days  # 平均
            self.cvd_df['Rt'] = power10_df(self.cvd_df['Rt_5_dB'])  # RtをdBから真値へ変換
            self.cvd_df['New_cases_rolling'] = power10_df(self.cvd_df['New_cases_rolling_dB'])  # 真値へ変換
        elif self.sw == 2:  # 実験的、未完
            self.cvd_df['New_cases_db'] = dB10_df(self.cvd_df['New_cases'])  # dB
            self.cvd_df['New_cases_db_diff1'] = self.cvd_df['New_cases'].diff(periods=7) / 7  # 1日差分
            self.cvd_df['New_cases_rolling'] = self.cvd_df['New_cases_db_diff1'].cumsum()
        #
        self.cvd_df['Rt'] = self.cvd_df['Rt'].replace([-np.Inf, np.Inf, 0.0], np.NaN)  # データクリーニング
        self.cvd_df['i14'] = self.cvd_df['New_cases'] * (self.cvd_df['Rt'] ** 2.8)  # 2.8=14/5
        self.cvd_df['i14'] = (np.floor((self.cvd_df['i14'].replace([-np.Inf, np.Inf, np.NaN], 0.0) + 0.5))).astype(
            'int')
        self.cvd_df['date'] = self.cvd_df.index  # 日付をいれる(前に動かしてはいけない)
        self.cvd_df.to_csv(self.dir + self.area + '.csv')  # CSVへ保存する
        # print('compute_rt\n')
        # print(self.cvd_df)
        return self

    # (5)感染者数の生データを分離してCSVへ保存-NHKデータオブジェクト
    def nhk_separete_data(self):  # for NHK
        # NHKデータは県別ではなく、一緒に入っている
        # NHKデータを読み、'都道府県名.csv'、'date'、'new_cases'を作る
        # 都道府県名を抜き出し、都道府県名.csv へ保存する
        print('nhk_separete_data')
        self.df_df['日付'] = self.df_df['日付'].replace('/', '-', regex=True)  # 2020/01/01 から2020-01-01　へ変換
        self.prefecture_names_se = self.df_df['都道府県名']  # 件名リストを作る
        self.prefecture_names_se = self.prefecture_names_se.drop_duplicates().reset_index(
            drop=True)  # 都道府県名の重複を削除する#inexを連番に振りなおす
        self.df_df.set_index('日付', drop=True)  # 日付をindexにセットする
        self.prefecture_names_se.to_csv(self.dir + '都道府県名.csv', index=False)  # 都道府県名をCSVファイルに保存する
        for prefecture_name in self.prefecture_names_se:  # 全県を回して、dateとNew_casesをcvd_dfへコピーする
            temp: object = self.df_df[self.df_df['都道府県名'] == prefecture_name]  # 都道府県名でデータを抜き出す
            self.cvd_df = pd.DataFrame(columns=self.col)  # データクリア
            self.cvd_df[['date', 'New_cases']] = temp[['日付', '各地の感染者数_1日ごとの発表数']].reset_index(
                drop=True)  # 日付、１日の感染者数をコピー、インデックスをリセット
            self.cvd_df['New_cases'].fillna(0, inplace=True)  # 新規感染者数のデータが無い所は０と置く
            self.cvd_df.to_csv('NHK全国データ/' + prefecture_name + '.csv', index=False)  # 県別にデータを保存
        return self

    # (6)感染者数の生データを県別compute_rt処理-NHKデータオブジェクト
    def nhk_compute_rt(self):  # for NHK
        print('nhk_compute_rt')
        for prefecture_name in self.prefecture_names_se:  # 全県を回す
            self.cvd_df = pd.DataFrame()  # クリア
            csv_file_name = self.dir + prefecture_name + '.csv'
            self.cvd_df = pd.read_csv(csv_file_name)  # 県別リストを読む
            self.cvd_df = self.cvd_df.set_index('date', drop=True)  # indexセット
            self.compute_rt()  # 実効再生産数を計算本体
            self.cvd_df.to_csv(csv_file_name, index=False)  # 県別リストを保存
        return self

    # (7)県別CSVから感染者数i,実効再生産数Rtエクセルへシートを分けて保存-NHKデータオブジェクト
    def nhk_make_excel_sheet(self):  # for NHK\
        print('nhk_make_excel_sheet')
        df_df = pd.read_csv('./' + self.file_name)  # df_dfを読む
        self.rt_all_df = pd.DataFrame(dtype=object)  # rt_all_df,全県のRt
        self.new_cases_all_df = pd.DataFrame(dtype=object)  # rt_all_df,全県のRt
        self.i14_all_df = pd.DataFrame(dtype=object)  # rt_all_df,全県のRt
        self.rt_all_df['date'] = df_df['日付'].copy  # 日付をコピー
        self.rt_all_df = self.rt_all_df.set_index('date', drop=True)  # indexセット
        for prefecture_name in self.prefecture_names_se:  # 全県を回す
            cvd_df: object = pd.read_csv(self.dir + prefecture_name + '.csv')  # 県別リストを保存
            cvd_df = cvd_df.set_index('date', drop=True)  # indexセット
            self.rt_all_df[prefecture_name] = cvd_df['Rt']  # Rtを保存
            self.new_cases_all_df[prefecture_name] = cvd_df['New_cases']  # Rtを保存
            self.i14_all_df[prefecture_name] = cvd_df['i14']  # Rtを保存
        # print('i14', self.i14_all_df)
        with pd.ExcelWriter('Rt.xlsx') as writer:  # 書き込み
            self.rt_all_df.to_excel(writer, engine='openpyxl', sheet_name='Rt')
            self.new_cases_all_df.to_excel(writer, engine='openpyxl', sheet_name='New_cases')
            self.i14_all_df.to_excel(writer, engine='openpyxl', sheet_name='i14')
        return self

    # (8)感染者数と実効再生産数のソートを追加-NHKデータオブジェクト
    def nhk_do_sort(self):  # for NHK
        print('nhk_do_sort')
        rt_df: object = pd.read_excel('Rt.xlsx', sheet_name='Rt', index_col='date')  # 県別Rtを読む
        rt_sorted_channel, rt_sorted_value = self.sort_by_values_and_return_columns(rt_df)  # Rtを県別にソートする
        new_cases_df: object = pd.read_excel('Rt.xlsx', sheet_name='New_cases')  # 県別Rtを読む
        # print('new_cases_df\n',new_cases_df)
        new_cases_df = new_cases_df.set_index('date', drop=True)  # indexを日付にセット
        new_cases_sorted_channel, new_cases_sorted_value = self.sort_by_values_and_return_columns(
            new_cases_df)  # Rtを県別にソートする
        new_cases_sorted_channel.index = new_cases_df.index  # indexを日付にセット
        new_cases_sorted_value.index = new_cases_df.index  # indexを日付にセット
        # print('new_cases_sorted_channel\n',new_cases_sorted_channel)
        # print('new_cases_sorted_value\n',new_cases_sorted_value)
        with pd.ExcelWriter('Rt.xlsx') as writer:  # 書き込み
            new_cases_sorted_channel.to_excel(writer, engine='openpyxl', sheet_name='ソートずみ新規感染者数県別順位')
            new_cases_sorted_value.to_excel(writer, engine='openpyxl', sheet_name='ソートずみ新規感染者数')
            rt_sorted_channel.to_excel(writer, engine='openpyxl', sheet_name='ソートずみ実効再生産数県別順位')
            rt_sorted_value.to_excel(writer, engine='openpyxl', sheet_name='ソートずみ実効再生産数')
            #
            self.new_cases_all_df.to_excel(writer, engine='openpyxl', sheet_name='新規感染者数')
            self.rt_all_df.to_excel(writer, engine='openpyxl', sheet_name='実効再生産数')
            self.rt_all_df.to_excel(writer, engine='openpyxl', sheet_name='Rt')
            self.new_cases_all_df.to_excel(writer, engine='openpyxl', sheet_name='New_cases')
            self.i14_all_df.to_excel(writer, engine='openpyxl', sheet_name='i14')
        return self

    # (9)ワクチン接種割合を追加-NHKデータオブジェクト
    def nhk_add_Ratio(self):  # for NHK-data
        print('nhk_add_Ratio')
        # self.vaccine_df = pd.read_excel('ワクチン接種回数2.xlsx')  # 2020-04-01形式、エクセル禁止、2020-4-1はダメ
        self.vaccine_df = pd.read_excel('ワクチン/ワクチン.xlsx')  # 2020-04-01形式、エクセル禁止、2020-4-1はダメ
        self.vaccine_df.index = pd.to_datetime(self.vaccine_df['日付'])  # CSVデータを読み込む
        self.vaccine_df['Ratio1'] = self.vaccine_df['累積1回目%']
        self.vaccine_df['Ratio2'] = self.vaccine_df['累積2回目%']
        # import os# 削除
        # try:
        #    os.removedirs('./NHK全国データ')
        #    os.mkdir('./NHK全国データ')
        # except IOError:
        #    pass
        for self.prefecture_name in self.prefecture_names_se:  # 県を回す
            print('nhk_add_Ratio:self.prefecture_name', self.prefecture_name)
            csv_name = self.dir + self.prefecture_name + '.csv'
            self.cvd_df = pd.read_csv(csv_name)  # 2020-04-01形式、エクセル禁止、2020-4-1はダメ
            self.cvd_df.index = pd.to_datetime(self.cvd_df['date'])  #
            self.cvd_df['Ratio1'] = self.vaccine_df['Ratio1']  # データを取り出す
            self.cvd_df['Ratio1'].fillna(method='ffill', inplace=True)
            self.cvd_df['Ratio2'] = self.vaccine_df['Ratio2']  # データを取り出す
            self.cvd_df['Ratio2'].fillna(method='ffill', inplace=True)
            self.cvd_df['累積1回目'] = self.vaccine_df['累積1回目']  # データを取り出す
            self.cvd_df['累積1回目'].fillna(method='ffill', inplace=True)
            self.cvd_df['累積2回目'] = self.vaccine_df['累積2回目']  # データを取り出す
            self.cvd_df['累積2回目'].fillna(method='ffill', inplace=True)
            self.cvd_df.to_csv(csv_name)

    # (10)NHKのグラフを書くオブジェクト# 全県を回す
    def nhk_plot_data(self):  # for NHK-data
        print('nhk_plot_data')
        new_cases_sorted_channel = pd.read_excel('Rt.xlsx', sheet_name='ソートずみ新規感染者数県別順位',
                                                 index_col='date', engine='openpyxl')
        self.emp_se = new_cases_sorted_channel.iloc[-1]  # 最後の行をとる
        new_cases_sorted_value = pd.read_excel('Rt.xlsx', sheet_name='ソートずみ新規感染者数', engine='openpyxl', index_col='date')
        self.empv_se = new_cases_sorted_value.iloc[-1]  # 最後の行をとる
        for self.emergency_prefectures_number, self.prefecture_name in enumerate(self.emp_se):  # 県を回す
            self.plot_data()  # プロット
        return self

    # (11)グラフを描くオブジェクト#1県だけ
    def plot_data(self):  # （1）、（2）グラフ

        # import matplotlib.pyplot as plt
        # import japanize_matplotlib

        import japanize_matplotlib
        self.cvd_df = pd.DataFrame(dtype=object)  # クリア
        self.cvd_df = pd.read_csv(self.dir + self.prefecture_name + '.csv')  # 2020-04-01形式、エクセル禁止、2020-4-1はダメ
        # fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        figure = plt.subplots(1, 1, figsize=(12, 8))
        x0 = pd.to_datetime(self.cvd_df['date'])  # CSVデータを読み込む
        y0 = self.cvd_df['New_cases']  # データを取り出す
        y1 = self.cvd_df['New_cases_rolling']  # データを取り出す
        y3 = self.cvd_df['Rt']  # データを取り出す
        # y4 = self.cvd_df['Ratio1']  # データを取り出す
        # y5 = self.cvd_df['Ratio2']  # データを取り出す
        plt.rcParams["font.size"] = 12
        if self.sw == 0:
            # plt.semilogy(x0, y0, '.r', x0, y1, '-b', x0, y3, '.-k', x0, y4, '.-m', x0, y5, '.-g',
            #             linewidth=3.0)  # 対数でグラフを描く
            plt.semilogy(x0, y0, '.r', x0, y1, '-b', x0, y3, '.-k',  # x0, y4, '.-m', x0, y5, '.-g',
                         linewidth=3.0)  # 対数でグラフを描く
            plt.ylim([1E-01, 1E04])  # ｙ軸の範囲
        elif self.sw != 0:
            plt.semilogy(x0, y0, '.r', x0, y1, '-b', x0, y3, '-k',  # x0, y4, '.-m', x0, y5, '.-g',
                         linewidth=2.5)  # 対数でグラフを描く
            plt.ylim([1E-01, 1E04])  # ｙ軸の範囲
        plt.xlabel("date:日付")  # 軸のラベル
        # plt.ylabel("i[人], iave(7)[人], Rt, rv1[%], rv2[%]")  # 軸のラベル
        plt.xlim(pd.to_datetime(["2020-06-01", "2021-10-01"]))  # ｘ軸の範囲

        #  新規感染者数New_cases,同7日平均データNew_cases_rolling,実効再生産数Rtを抜き出して、文字列へ変換する
        date_list = self.cvd_df['date'].to_list()  # 日付をデータフレームからリストへ変換
        self.date_str = "{}".format(date_list[-1])  # 日付を文字列へ変換
        t0_list = self.cvd_df['New_cases'].to_list()  # 新規感染者数をデータフレームからリストへ変換
        t0 = "{:,.0f}".format(t0_list[-1])  # 最後の新規感染者数文字列へ変換
        t1_list = self.cvd_df['New_cases_rolling'].to_list()  # 7日平均をデータフレームからリストへ変換
        t1 = "{:.1f}".format(t1_list[-1])  # 最後の7日平均を文字列へ変換

        t3_list = self.cvd_df['Rt'].to_list()  # 実効再生産数をデータフレームからリストへ変換
        t3 = "{:.2f}".format(t3_list[-1])  # 最後の実効再生産数を文字列へ変換

        t4_list = self.cvd_df['Ratio1'].to_list()  # ワクチン接種率をデータフレームからリストへ変換
        t4 = "{:.2f}".format(t4_list[-1])  # 最後のワクチン接種率を文字列へ変換

        t5_list = self.cvd_df['Ratio2'].to_list()  # ワクチン接種率をデータフレームからリストへ変換
        t5 = "{:,.2f}".format(t5_list[-1])  # 最後のワクチン接種率を文字列へ変換

        t6_list = (self.cvd_df['累積1回目'] / 10000).to_list()  # ワクチン接種回数pをデータフレームからリストへ変換
        t6 = "{:,.0f}".format(t6_list[-1])  # 最後のワクチン接種回数を文字列へ変換

        t7_list = (self.cvd_df['累積2回目'] / 10000).to_list()  # ワクチン接種回数をデータフレームからリストへ変換
        t7 = "{:,.0f}".format(t7_list[-1])  # 最後のワクチン接種回数を文字列へ変換

        plt.legend(['i:新規感染者数(' + t0 + '人)',
                    'iave(7):同7日平均(' + t1 + '人)',
                    'Rt:(日毎)実効再生産数(' + t3 + ')', ])  # 数値付き凡例を表示
        #            'rv1:全国ワクチン1回接種' + t6 + '万回(' + t4 + '%)',
        #            'rv2:全国ワクチン2回接種' + t7 + '万回(' + t5 + '%)']

        plt.grid(ls='--', which="both")  # 罫線

        # ファイルに保存
        plt.title(self.prefecture_name + '    ' + self.date_str)  # タイトル
        plt.savefig(self.dir + str(self.emergency_prefectures_number) + self.prefecture_name + '.png')
        plt.show()  # グラフを表示
        return self

    # (12)ワクチン接種の生データをダウンロード-首相官邸★
    # class Vaccine(Covid19):
    def vaccine_iltupan(self):  # 一般
        temp = pd.read_excel('./vaccination_data5.xlsx', sheet_name='一般接種',
                             skiprows=[0, 1, 3, 4, 5], header=0, dtype=str)  # データを読む
        iltupan_df = pd.DataFrame()
        iltupan_df['日付'] = temp.iloc[:, 0]
        iltupan_df['曜日'] = temp.iloc[:, 1]
        iltupan_df['一般1回目_ファイザー'] = temp.iloc[:, 3].astype('float')
        iltupan_df['一般1回目_武田モデルナ'] = temp.iloc[:, 4].astype('float')
        iltupan_df['一般1回目_アストラゼネカ'] = temp.iloc[:, 5].astype('float')
        #
        iltupan_df['一般2回目_ファイザー'] = temp.iloc[:, 6].astype('float')
        iltupan_df['一般2回目_武田モデルナ'] = temp.iloc[:, 7].astype('float')
        iltupan_df['一般2回目_アストラゼネカ'] = temp.iloc[:, 8].astype('float')
        #
        iltupan_df['高齢1回目_ファイザー'] = temp.iloc[:, 9].astype('float')
        iltupan_df['高齢1回目_武田モデルナ'] = temp.iloc[:, 10].astype('float')
        iltupan_df['高齢1回目_アストラゼネカ'] = temp.iloc[:, 11].astype('float')
        #
        iltupan_df['高齢2回目_ファイザー'] = temp.iloc[:, 12].astype('float')
        iltupan_df['高齢2回目_武田モデルナ'] = temp.iloc[:, 13].astype('float')
        iltupan_df['高齢1回目_アストラゼネカ'] = temp.iloc[:, 14].astype('float')
        #
        iltupan_df['一般1回目'] = iltupan_df['一般1回目_ファイザー']\
            .add(iltupan_df['一般1回目_武田モデルナ'], fill_value=0)\
            .add(iltupan_df['一般1回目_アストラゼネカ'], fill_value=0)
        #
        iltupan_df['一般2回目'] = iltupan_df['一般2回目_ファイザー']\
            .add(iltupan_df['一般2回目_武田モデルナ'], fill_value=0)\
            .add(iltupan_df['一般2回目_アストラゼネカ'], fill_value=0)
        #
        iltupan_df['一般接種回数'] = iltupan_df.iloc[:, 2]  # 一般接種回数をコピー
        iltupan_df['高齢1回目'] = iltupan_df['高齢1回目_ファイザー'].add(iltupan_df['高齢1回目_武田モデルナ'], fill_value=0)
        iltupan_df['高齢2回目'] = iltupan_df['高齢2回目_ファイザー'].add(iltupan_df['高齢2回目_武田モデルナ'], fill_value=0)
        iltupan_df['65歳未満1回目'] = iltupan_df['一般1回目'].sub(iltupan_df['高齢1回目'], fill_value=0)
        iltupan_df['65歳未満2回目'] = iltupan_df['一般2回目'].sub(iltupan_df['高齢2回目'], fill_value=0)
        iltupan_df.dropna(axis='index', subset=['曜日'], inplace=True)  # NaNのある行を削除
        iltupan_df.to_csv('./ワクチン/0一般・高齢0.csv')  # 保存
        self.vc_df = iltupan_df
        self.vc_df.to_csv('./ワクチン/0一般・高齢1.csv', index=False)  # 保存

    # 医療ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
    def vaccine_iryou(self):
        temp = pd.read_excel('./vaccination_data5.xlsx',
                             sheet_name='医療従事者等', skiprows=[0, 1, 3, 4], dtype=str)
        temp.dropna(axis='columns', how='all', inplace=True)  # NaNのある列を削除
        temp.dropna(axis='index', subset=['曜日'], inplace=True)  # NaNのある行を削除
        iryou_df = pd.DataFrame()
        iryou_df['日付'] = temp.iloc[:, 0]
        iryou_df['曜日'] = temp.iloc[:, 1]
        iryou_df['医療1回目_ファイザー'] = temp.iloc[:, 3].astype('float')
        iryou_df['医療1回目_武田モデルナ'] = temp.iloc[:, 4].astype('float')
        iryou_df['医療2回目_ファイザー'] = temp.iloc[:, 5].astype('float')
        iryou_df['医療2回目_武田モデルナ'] = temp.iloc[:, 6].astype('float')
        iryou_df['医療1回目'] = iryou_df['医療1回目_ファイザー'].add(iryou_df['医療1回目_武田モデルナ'], fill_value=0)
        iryou_df['医療2回目'] = iryou_df['医療2回目_ファイザー'].add(iryou_df['医療2回目_武田モデルナ'], fill_value=0)
        iryou_df['医療接種回数'] = temp.iloc[:, 2].astype('float')  # 医療接種回数をコピー

        iryou_df.to_csv('./ワクチン/1医療.csv', index=False,encoding='s-jis')  # 保存
        self.vc_df = self.vc_df.merge(iryou_df, how='outer', on='日付')

    # 職域ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
    def vaccine_syokuiki(self):  # 職域
        temp_df = pd.read_excel('./vaccination_data5.xlsx',
                                sheet_name='職域接種',header=0, skiprows=[0, 1, 3], dtype=str)  # データを読む
        temp_df.dropna(axis='columns', how='all', inplace=True)  # NaNのある列を削除
        temp_df.dropna(axis='index', subset=['曜日'], inplace=True)  # NaNのある行を削除
        #temp_df.to_csv('./ワクチン/1職域0.csv',encoding='s-jis')  # 保存
        temp2=pd.DataFrame()
        temp2['日付'] = temp_df.iloc[:, 2]
        temp2['職域1回目_累積'] = temp_df.iloc[:, 5].astype('float')
        temp2['職域2回目_累積'] = temp_df.iloc[:, 6].astype('float')
        temp2['職域接種回数_累積'] = temp_df.iloc[:, 4].astype('float')
        #temp_df.to_csv('./ワクチン/1職域1.csv',encoding='s-jis')  # 保存
        syokuiki_df = pd.DataFrame()
        syokuiki_df['日付'] = self.vc_df['日付']
        syokuiki_df = syokuiki_df.merge(temp2, how='outer', on='日付')
        syokuiki_df = syokuiki_df.sort_index(ascending=False)
        syokuiki_df['職域1回目_累積'].ffill(inplace=True)
        syokuiki_df['職域2回目_累積'].ffill(inplace=True)
        syokuiki_df['職域接種回数_累積'].ffill(inplace=True)
        syokuiki_df.to_csv('./ワクチン/1職域2.csv',encoding='s-jis')  # 保存
        syokuiki_df['職域1回目'] = syokuiki_df['職域1回目_累積'].fillna(0).diff()
        syokuiki_df['職域2回目'] = syokuiki_df['職域2回目_累積'].fillna(0).diff()
        syokuiki_df['職域接種回数'] = syokuiki_df['職域接種回数_累積'].fillna(0).diff()
        syokuiki_df = syokuiki_df.sort_index(ascending=True)
        syokuiki_df.to_csv('./ワクチン/1職域3.csv')  # 保存
        self.vc_df = self.vc_df.merge(syokuiki_df, how='outer', on='日付')
        self.vc_df.to_csv('./ワクチン/1職域4.csv', index=False,encoding='s-jis')  # 保存

    # 累積ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
    def vaccine_ruiseki(self):
        self.vc_df['65歳未満1回目'] = self.vc_df['65歳未満1回目'].add(self.vc_df['職域1回目'], fill_value=0)
        self.vc_df['65歳未満2回目'] = self.vc_df['65歳未満2回目'].add(self.vc_df['職域2回目'], fill_value=0)
        self.vc_df['接種回数'] = self.vc_df['一般接種回数'].add(self.vc_df['医療接種回数'], fill_value=0).add(self.vc_df['職域接種回数'],
                                                                                              fill_value=0)  # 加算
        self.vc_df['全体1回目'] = self.vc_df['一般1回目'].add(self.vc_df['医療1回目'], fill_value=0).add(self.vc_df['職域1回目'],
                                                                                             fill_value=0)
        self.vc_df['全体2回目'] = self.vc_df['一般2回目'].add(self.vc_df['医療2回目'], fill_value=0).add(self.vc_df['職域2回目'],
                                                                                             fill_value=0)
        self.vc_df['1日接種回数'] = self.vc_df['全体1回目'].add(self.vc_df['全体2回目'], fill_value=0)
        """
        self.vc_df['65歳未満1回目'] = self.vc_df['65歳未満1回目']+self.vc_df['職域1回目']
        self.vc_df['65歳未満2回目'] = self.vc_df['65歳未満2回目']+self.vc_df['職域2回目']
        self.vc_df['接種回数'] = self.vc_df['一般接種回数'] + self.vc_df['医療接種回数'] + self.vc_df['職域接種回数']
        self.vc_df['全体1回目'] = self.vc_df['一般1回目'] + self.vc_df['医療1回目'] + self.vc_df['職域1回目']
        self.vc_df['全体2回目'] = self.vc_df['一般2回目'] + self.vc_df['医療2回目'] + self.vc_df['職域2回目']
        self.vc_df['1日接種回数'] = self.vc_df['全体1回目'] + self.vc_df['全体2回目']
        """
        # 　累積ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
        self.vc_df.sort_index(inplace=True, ascending=False)  # インデックスでソート
        self.vc_df['累積医療接種回数'] = self.vc_df['医療接種回数'].cumsum()  # 累積回数を積算
        self.vc_df['累積職域接種回数'] = self.vc_df['職域接種回数'].cumsum()  # 累積回数を積算
        self.vc_df['累積65歳未満1回目'] = self.vc_df['65歳未満1回目'].cumsum()  # 累積回数を積算
        self.vc_df['累積65歳未満2回目'] = self.vc_df['65歳未満2回目'].cumsum()  # 累積回数を積算
        self.vc_df['累積1回目'] = self.vc_df['全体1回目'].cumsum()  # 累積回数を積算
        self.vc_df['累積2回目'] = self.vc_df['全体2回目'].cumsum()  # 累積回数を積算
        self.vc_df['累積(全体)'] = self.vc_df['1日接種回数'].cumsum()  # 累積回数を積算を積算
        self.vc_df=self.vc_df.replace(np.nan,0)
        self.vc_df=self.vc_df.replace(np.inf,0)
        #

    # 重複ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
    def vaccine_cyoufuku(self):
        # 　「総接種回数シート」の「重複」を「累計(全体)」から減算するーーーーーーーーーー
        sou_df = pd.read_excel('./vaccination_data5.xlsx', sheet_name='総接種回数', skiprows=[0])  # データを読む
        sou_df.dropna(axis='columns', how='all', inplace=True)  # NaNのある列を削除
        sou_df['日付'] = sou_df['公表日']  # 集計日をindexにセット
        sou_df.set_index('日付', drop=True)  # 集計日をindexにセット
        sou_df = sou_df.replace('―', None)
        sou_df['重複'] = sou_df['重複'].astype('float')
        sou_df['職域接種'] = sou_df['職域接種'].astype('float')
        sou_df.to_csv('./ワクチン/4総接.csv',encoding='s-jis')  # 保存
        #
        self.vc_df['累積(全体)'] = self.vc_df['累積(全体)'].sub(sou_df['重複'], fill_value=0)
        ##
        # self.vc_df['累積(全体)'] +=1101698#補正
        ##

    # パーセント計算ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
    def vaccine_percent(self):
        self.vc_df['累積1回目%'] = self.vc_df['累積1回目'] * 100 / self.jinkou  # データを取り出す # 累積回数を積算
        self.vc_df['累積2回目%'] = self.vc_df['累積2回目'] * 100 / self.jinkou  # データを取り出す  # 累積回数を積算
        self.vc_df['累積(全体)%'] = self.vc_df[
                                    '累積(全体)'] * 100 / self.jinkou  # データを取り出す)  # 累積回数データを取り出す)  # 累積回数を積算
        self.vc_df['累積65歳未満1回目%'] = self.vc_df['累積65歳未満1回目'] * 100 / self.jinkou  # データを取り出す # 累積回数を積算
        self.vc_df['累積65歳未満2回目%'] = self.vc_df['累積65歳未満2回目'] * 100 / self.jinkou  # データを取り出す  # 累積回数を積算
        self.vc_df['実効再生産数_ワクチン効果'] = (100 - (self.vc_df['累積1回目%'] - self.vc_df['累積2回目%']) / 3 - self.vc_df[
            '累積2回目%'] * 0.94) / 100  #
        # self.vc_df.dropna(axis='index', subset=['累積(全体)'], inplace=True)  # NaNのある行を削除
        print('self.vc_df\n', self.vc_df.head())
        self.vc_df.to_excel('./ワクチン/ワクチン.xlsx', index=None)  # 保存

    # グラフを描くーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
    def vaccine_graph(self):
        import matplotlib.pyplot as plt
        import japanize_matplotlib
        fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(12, 6), dpi=100)  # フィガーをセット
        x0 = pd.to_datetime(self.vc_df['日付'])  # CSVデータを読み込む
        last_d = x0.iloc[-1].strftime("20%y-%m-%d")  # 最新の日付
        self.vc_df['1回目_1日接種回数_7日平均'] = self.vc_df['累積1回目'].diff(1)
        rest_days_to_finish = int(
            (self.jinkou * 0.50 - self.vc_df['累積1回目'].iloc[-1]) / self.vc_df['1回目_1日接種回数_7日平均'].iloc[
                -1])  # 1日接種回数_7日平均,2日前
        print('rest_days_to_finish', rest_days_to_finish)
        finish_d = (x0.iloc[-1] + pd.offsets.Day(rest_days_to_finish)).strftime("%-m月%-d日")
        # ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
        y6 = self.vc_df['累積(全体)%']  # データを取り出す
        y7 = self.vc_df['累積1回目%']  # データを取り出す
        y8 = self.vc_df['累積2回目%']  # データを取り出す
        y9 = self.vc_df['累積65歳未満1回目%']  # データを取り出す
        y10 = self.vc_df['累積65歳未満2回目%']  # データを取り出す
        #
        v_7days_average = "{:,.0f}".format(self.vc_df['1回目_1日接種回数_7日平均'].iloc[
                                               -3])
        v_all = "{:,.0f}".format(self.vc_df['累積(全体)'].iloc[-1])
        v_1_dose = "{:,.0f}".format(self.vc_df['累積1回目'].iloc[-1])
        v_2_dose = "{:,.0f}".format(self.vc_df['累積2回目'].iloc[-1])
        v_1_percent = "{:.1f}".format(self.vc_df['累積1回目%'].iloc[-1])
        v_2_percent = "{:.1f}".format(self.vc_df['累積2回目%'].iloc[-1])
        v_l65_1_percent = "{:.1f}".format(self.vc_df['累積65歳未満1回目%'].iloc[-1])
        v_l65_2_percent = "{:.1f}".format(self.vc_df['累積65歳未満2回目%'].iloc[-1])

        v_jinkou = "{:,.0f}".format(self.jinkou)
        # plt2ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
        ax2.set_xlim(pd.to_datetime(["2021-06-1", "2021-10-1"]))  # ｘ軸の範囲
        ax2.set_ylim([0.1, 100])  # ｙ軸の範囲
        ax2.grid(ls='--', which="both")  # 罫線
        ax2.set_xlabel("date:日付", fontsize=10)  # 軸のラベル
        ax2.set_ylabel("人口に対する接種割合[％]", fontsize=10)  # 軸のラベル
        ax2.text(pd.to_datetime(['2021-06-05']), 30, '黒線：接種割合50%\n人口' + v_jinkou + '人(令和2年1月1日,外国人含む)',
                 color='black')  # バーグラフを描く
        # ax2.semilogy(x0, y7, color='red', label="1回目%(" + v_1_percent + ")")  # グラフを描く
        # ax2.semilogy(x0, y8, color='blue', label="2回目%(" + v_2_percent + ")")  # グラフを描く
        # ax2.semilogy(pd.to_datetime(["2021-04-1", "2021-08-1"]), [50, 50], color='black', linewidth=0.5)  # グラフを描く
        # ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
        ax2.plot(x0, y7, color='red', label="1回目%(" + v_1_percent + ")")  # グラフを描く
        ax2.plot(x0, y8, color='blue', label="2回目%(" + v_2_percent + ")")  # グラフを描く
        ax2.plot(x0, y9, color='cyan', label="1回目%(" + v_l65_1_percent + ")")  # グラフを描く
        ax2.plot(x0, y10, color='magenta', label="2回目%(" + v_l65_2_percent + ")")  # グラフを描く
        ax2.plot(pd.to_datetime(["2021-06-1", "2021-10-1"]), [50, 50], color='black', linewidth=0.5)  # グラフを描く
        # plt1ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
        y1 = self.vc_df['全体1回目'] / 10000  # データを取り出す
        y2 = self.vc_df['全体2回目'] / 10000  # データを取り出す
        ax1.set_xlim(pd.to_datetime(["2021-06-1", "2021-10-1"]))  # ｘ軸の範囲
        ax1.set_ylim([0, 1000])  # ｙ軸の範囲
        ax1.set_ylabel("接種回数[万回]", fontsize=10)  # 軸のラベル
        # 文字設定ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
        text_header = '新型コロナウイルス･ワクチン接種回数'
        text_sou_kaisuu = "総接種回数 " + v_all + '回'
        text_ichi_nichi_kaisuu = "1回接種回数7日平均" + v_7days_average + '回'
        text_1_ruikei = "1回接種済(" + v_1_dose + '人'
        text_1_percent = ',人口比' + v_1_percent + "%)"
        text_2_ruikei = "2回接種済(" + v_2_dose + '人'
        text_2_percent = ',人口比' + v_2_percent + "%)"
        text_jinkou = '日本の人口(令和2年1月1日,外国人含む)' + v_jinkou + '人'
        text_days = '1回接種済人口比50％予想日' + finish_d
        text_footer = last_d + '閲覧,データ:首相官邸HP'
        text_all = pd.Series(text_header + '\n\n' +
                             text_sou_kaisuu + '\n' +
                             text_ichi_nichi_kaisuu + '\n' +
                             text_1_ruikei + text_1_percent + '\n' +
                             text_2_ruikei + text_2_percent + '\n' +
                             text_jinkou + '\n' +
                             text_days + '\n' +
                             text_footer)
        print(text_all)
        text_all.to_csv('./ワクチン/ワクチン接種.txt', index=False, header=False, sep=' ', quotechar=' ')
        # グラフを描く
        ax1.bar(x0, y2, color='blue', bottom=y1, label=text_2_ruikei + text_2_percent)  # 2回目
        ax1.bar(x0, y1, color='red', label=text_1_ruikei + text_1_percent)  # 1回目
        ax1.text(pd.to_datetime(['2021-06-05']), 80, text_sou_kaisuu + '\n' + text_ichi_nichi_kaisuu,
                 color='black')  # 総接種回数
        ax1.set_title(text_header + '  ' + text_footer, fontsize=10)
        ax1.grid(ls='--', which="both")  # 罫線
        ax1.legend()
        fig.savefig('./ワクチン/ワクチン接種.png')
        fig.show()

    def kantei_vaccine_dose_data(self):  #
        print('kantei_download_data')
        self.vaccine_iltupan()  # 一般
        self.vaccine_iryou()  # 医療
        self.vaccine_syokuiki()  # 職域
        self.vaccine_ruiseki()  # 累積
        self.vaccine_cyoufuku()  # 重複
        self.vaccine_percent()  # パーセント計算
        self.vaccine_graph()  # グラフ
        return self

    # (13)データフレームdf_dfから、もっとも大きいデータとコードを取り出す関数
    @staticmethod
    def sort_by_values_and_return_columns(df_df):
        # １行ごとにデータを取り出し（日付順）
        # 列方向に値でソートし
        # 値でソートされたカラム名（項目1、項目2、項目3）を返す
        #
        # 　（１）元データdf_dfの構造
        #  日付　　　　   項目1　項目2　項目3　....
        #  2019/06/02　　1.02　　2.01　　1.04　　value
        #  2019/06/01　　1.32　　2.41　　1.43　

        # 　（2）転置
        #  日付   2019/06/02　2019/06/01　　　....
        #  項目1 1.02　　      1.32　　　　　　
        #  項目2 2.01　　      2.41　
        #  項目3 1.04　　      1.43

        # 　（3）1列取り出す（この場合は１列目をSeriesへコピー）
        #  日付　　　　2019/06/02
        #  項目1　　1.02　　　　　　
        #  項目2　　2.01　
        #  項目3   1.04

        # 　（4）値でソート（１列のみ）
        #  日付　　　2019/06/02
        #  項目2　　2.01　
        #  項目3   1.04
        #  項目1　　1.02　　　　　　

        # 　（5）整列した項目を返す（1列のみ）、この処理を全ての列で繰り返す
        #  順番　　　　2019/06/02
        #  ０        項目2　
        #  1        項目3
        #  2        項目1　

        # 　（6）全ての日付でソートを行った結果をDataFrameへコピー
        #  順番　　　　2019/06/02　2019/06/01　　　....
        #  ０        項目2    項目2
        #  1        項目3    項目3
        #  2        項目1　　 項目1　　　　

        # 　（7）転置して日付毎の、もっとも強い項目順の並べかえ完成。resへコピー、リターン
        #  順番           0   1   2
        #  2019/06/02　項目2　項目3 　項目1
        #  2019/06/01　項目2　項目3 　項目1
        #
        # df_df = pd.DataFrame([[1, 3, 7], [4, 8, 5],
        #                      [9, 5, 3], [6, 9, 10], [11, 2, 1]])
        # df_df.index = ["1/1", "1/2", "1/3", "1/14", "1/16"]
        # df_df.columns = ["F1", "F2", "F3"]
        # print(df_df.T)
        # df_df.T
        #    1/1    1/2  1/3  1/14  1/16
        # F1    1    4    9     6    11
        # F2    3    8    5     9     2
        # F3    7    5    3    10     1
        # sorted_value
        #   1/1  1/2  1/3  1/14  1/16
        # 0    7    8    9    10    11
        # 1    3    5    5     9     2
        # 2    1    4    3     6     1
        # sorted_channel
        #  1/1 1/2 1/3 1/14 1/16
        # 0  F3  F2  F1   F3   F1
        # 1  F2  F3  F2   F2   F2
        # 2  F1  F1  F3   F1   F3

        df_df = df_df.T  # 転置する
        sorted_channel = pd.DataFrame()  # 結果保存用
        sorted_value = pd.DataFrame()  # 結果保存用
        for day in df_df.columns:  # （転置しているので）インデックス＝日付を全て回す
            df_ser = df_df[day].copy()
            temp = df_ser.sort_values(ascending=False, inplace=None)  # 1行をソートしカラムをとりだして結果を保存
            # print(day, temp)
            sorted_channel[day] = temp.index
            sorted_value[day] = temp.values
        return sorted_channel.T, sorted_value.T

    # (15)感染者数のカーブフィットするオブジェクト2
    def curve_fit_rt(self):  # Rtを最小計二乗法で計算する
        self.prefecture_name = '東京'
        self.cvd_df = pd.read_csv(self.dir + self.prefecture_name + '.csv')  # 2020-04-01形式、エクセル禁止、2020-4-1はダメ
        self.cvd_df['New_cases_dB'] = dB10_df(self.cvd_df['New_cases'])  # 区間毎に2次関数で区間毎に最小自乗近似
        self.cvd_df['New_cases_max'] = max_all_range(self.cvd_df['New_cases_dB'], t0=7)  # dB
        self.cvd_df['New_cases_median'] = median_all_range(self.cvd_df['New_cases_dB'], t0=7)  # dB
        self.cvd_df['New_cases_min'] = min_all_range(self.cvd_df['New_cases_dB'], t0=7)  # dB
        self.cvd_df['New_cases_rms2'] = fit_by_linear_all_range(self.cvd_df['No'], self.cvd_df['New_cases_median'],
                                                                t_segment=7, t_extension=2)  # 区間毎に2次関数で区間毎に最小自乗近似
        self.cvd_df['Rt7_dB'] = (self.cvd_df['New_cases_rms2'].diff(periods=7)) * (5 / 7)  # dB領域で差分をとりdB(実効再生産数Rt)を計算
        import matplotlib.pyplot as plt
        import japanize_matplotlib
        self.cvd_df.fillna(0, inplace=True)
        plt.plot(self.cvd_df['No'], self.cvd_df['Rt7_dB'], '-m',
                 self.cvd_df['No'], self.cvd_df['New_cases_dB'], '.r',
                 self.cvd_df['No'], self.cvd_df['New_cases_max'], '-k',
                 self.cvd_df['No'], self.cvd_df['New_cases_median'], '-k',
                 self.cvd_df['No'], self.cvd_df['New_cases_min'], '-k',
                 self.cvd_df['No'], self.cvd_df['New_cases_rms2'], '-k',
                 linewidth=0.5
                 # self.cvd_df['No'], self.cvd_df['New_cases_rms'], '-r'
                 )
        plt.grid(True)
        plt.savefig(self.dir + str(self.prefecture_name) + '_curve_fit.png')
        plt.show()
        self.cvd_df.to_csv(self.dir + str(self.prefecture_name) + '_curve_fit.csv')
        return self

    # (16-0)OWIDのグループバイオブジェクトを読む
    def our_world_in_data_goupby_read(self):
        # 新型コロナウイルス
        # owid-covid-data.xlsxをダウンロードし加工してプロットする
        # https://github.com/owid/covid-19-data/blob/master/public/data/owid-covid-data.xlsx
        import matplotlib.pyplot as plt
        print('owid_covid_data')
        tmp_df = pd.read_excel('./ワクチン/owid-covid-data.xlsx')
        tmp_df2 = tmp_df[['date',
                          'location',
                          'new_cases',
                          'new_cases_smoothed',
                          'new_cases_smoothed_per_million',
                          'reproduction_rate',
                          'people_vaccinated',
                          'people_fully_vaccinated',
                          'new_vaccinations',
                          'new_vaccinations_smoothed',
                          'total_vaccinations_per_hundred',
                          'people_vaccinated_per_hundred',
                          'people_fully_vaccinated_per_hundred',
                          'new_vaccinations_smoothed_per_million']]
        tmp_df2.to_pickle('./ワクチン/owid.pkl')  # 保存
        return self

    # (16)OWIDのグループバイオブジェクト
    def covid_goupby(self):
        import matplotlib.pyplot as plt
        print('owid_covid_data')
        # self.our_world_in_data_goupby_read()
        # return
        vc_df = pd.read_pickle('./ワクチン/owid.pkl')  # 保存
        group_by = vc_df.groupby('location')
        x0 = pd.to_datetime(group_by.get_group('Japan')['date'])  # 日付に変換
        y01 = group_by.get_group('Japan')['new_cases_smoothed_per_million'] / 10  # 新規感染者数の７日平均
        y02 = group_by.get_group('Japan')['people_vaccinated_per_hundred']  # 累積１回ワクチン接種者数
        x1 = pd.to_datetime(group_by.get_group('Israel')['date'])  # 日付に変換
        y11 = group_by.get_group('Israel')['new_cases_smoothed_per_million'] / 10  # 新規感染者数の７日平均
        y12 = group_by.get_group('Israel')['people_vaccinated_per_hundred']  # 累積１回ワクチン接種者数
        plt.plot(x0, y01, '-g', x0, y02, '-m', x1, y11, '-b', x1, y12, '-r')
        plt.grid(ls='--', which="both")  # 罫線
        plt.savefig('./ワクチン/group_by.png')
        plt.show()
        return self

    # (17)Rt-i予測関数（未使用）
    def predict_i(self):
        def compute_rt2(t1, t2, i1, i2, u=5.0):
            rt = np.NaN
            if (np.abs(t2 - t1) > 0) and (i1 > 0) and (i2 > 0) and (u > 0):
                rt = 10 ** ((np.log10(i2) - np.log10(i1)) * (u / (t2 - t1)))
            return rt

        def compute_rt(t0, t, i0, i, u=5.0):
            return (i / i0) ** (u / (t - t0))
            # def compute_rt(t1, t2, i1, i2, u=5.0):
            #    rt = np.NaN
            #    if (np.abs(t2 - t1) > 0) and (i1 > 0) and (i2 > 0) and (u > 0):
            #        rt = 10 ** ((np.log10(i2) - np.log10(i1)) * (u / (t2 - t1)))
            #    return rt
            #
            #    def compute_i(t1, t2, i0, r0, u=5.0):
            #    return i0 * r0 ** ((t2 - t1) / u)

        def compute_i(t0, t, i0, rt, u=5.0):
            return i0 * rt ** ((t - t0) / u)

        def compute_i_mag(t, rt):
            u = 5
            return rt ** (t / u)

        t0, t, i0, rt = 0, 14, 1099, 1.34
        i3 = compute_i(t0, t, i0, rt)
        print('t0=', t0, 't=', t, ',i0=', i0, 'rt=', rt, 'i3=', round(i3))

        t0, t, i0, rt = 0, 14, 510, 1.17
        i3 = compute_i(t0, t, i0, rt)
        print('t0=', t0, 't=', t, ',i0=', i0, 'rt=', rt, 'i3=', round(i3))

        t0, t, i0, rt = 0, 14, 391, 1.36
        i3 = compute_i(t0, t, i0, rt)
        print('t0=', t0, 't=', t, ',i0=', i0, 'rt=', rt, 'i3=', round(i3))

        t0, t, i0, rt = 0, 14, 168, 1.38
        i3 = compute_i(t0, t, i0, rt)
        print('t0=', t0, 't=', t, ',i0=', i0, 'rt=', rt, 'i3=', round(i3))

    # (18)感染者数の生データを県別compute_i処理-NHKデータオブジェクト
    def nhk_compute_i(self):  # for NHK
        print('enter:nhk_compute_i')
        for prefecture_name in self.prefectures_se:  # 全県を回す
            # print(prefecture_name)
            self.cvd_df = pd.DataFrame()  # クリア
            csv_file_name = self.dir + prefecture_name + '.csv'
            self.cvd_df = pd.read_csv(csv_file_name)  # 県別リストを読む
            t, u = 14.0, 5.0  # t日後、世代間隔5日
            self.cvd_df['Mag14'] = self.cvd_df['Rt'] ** (t / u)  # t日後の新規感染者数を予想
            # self.cvd_df['i14'] = (self.cvd_df['New_cases'] * self.cvd_df['Mag14']).map(lambda x: float(decimal.Decimal(str(x)).
            #                decimal.quantize(Decimal('0'), rounding=ROUND_HALF_UP)))
            # pd.options.display.precision=0
            self.cvd_df['i14'] = pd.floor((self.cvd_df['New_cases'] * self.cvd_df['Mag14']) + 0.5)
            self.cvd_df.to_csv(csv_file_name, index=False)  # 県別リストを保存
        return self

    # (19)Twitter用文言CSVへ保存
    # '新規感染者数上位1-4位(i,Rt)\n#{}({}, {})\n#{}({}, {})\n#{}({}, {})\n#{}({}, {}, {})\n\n'.format(
    def make_text_for_twitter(self):
        print('make_text_for_twitter')
        new_cases_sorted_channel = pd.read_excel('Rt.xlsx', sheet_name='ソートずみ新規感染者数県別順位',
                                                 index_col='date', engine='openpyxl')
        self.emp_se = new_cases_sorted_channel.iloc[-1]  # 最後の行をとる
        new_cases_sorted_value = pd.read_excel('Rt.xlsx', sheet_name='ソートずみ新規感染者数', engine='openpyxl', index_col='date')
        self.empv_se = new_cases_sorted_value.iloc[-1]  # 最後の行をとる
        Rt = pd.read_excel('Rt.xlsx', sheet_name='Rt', engine='openpyxl', index_col='date')  # 県名対Rtデータベース
        Rt_se = Rt.iloc[-1]
        i14 = pd.read_excel('Rt.xlsx', sheet_name='i14', engine='openpyxl', index_col='date')  # 県名対Rtデータベース
        self.i14_se = i14.iloc[-1]

        # 　文言生成
        self.header1 = '新型コロナウイルス,新規感染者数:i0 1-4位\n2週間後の新規感染者数(予想):i14\n\n'
        self.header2 = '新型コロナウイルス,新規感染者数:i0 5-8位\n2週間後の新規感染者数(予想):i14\n\n'
        self.header3 = '新型コロナウイルス,新規感染者数:i0 9-12位\n2週間後の新規感染者数(予想):i14\n\n'
        self.header4 = '新型コロナウイルス,新規感染者数:i0 13-16位\n2週間後の新規感染者数(予想):i14\n\n'
        caption = "地域,      i0,  Rt,  i14\n"
        b = f"#{self.emp_se[0]:<s}, {self.empv_se[0]:>d}, {Rt_se[self.emp_se[0]]:>4.3f}, {self.i14_se[self.emp_se[0]]:>d}\n"
        c = f"#{self.emp_se[1]:<s}, {self.empv_se[1]:>d}, {Rt_se[self.emp_se[1]]:>4.3f}, {self.i14_se[self.emp_se[1]]:>d}\n"
        d = f"#{self.emp_se[2]:<s}, {self.empv_se[2]:>d}, {Rt_se[self.emp_se[2]]:>4.3f}, {self.i14_se[self.emp_se[2]]:>d}\n"
        e = f"#{self.emp_se[3]:<s}, {self.empv_se[3]:>d}, {Rt_se[self.emp_se[3]]:>4.3f}, {self.i14_se[self.emp_se[3]]:>d}\n"
        f = f"#{self.emp_se[4]:<s}, {self.empv_se[4]:>d}, {Rt_se[self.emp_se[4]]:>4.3f}, {self.i14_se[self.emp_se[4]]:>d}\n"
        g = f"#{self.emp_se[5]:<s}, {self.empv_se[5]:>d}, {Rt_se[self.emp_se[5]]:>4.3f}, {self.i14_se[self.emp_se[5]]:>d}\n"
        h = f"#{self.emp_se[6]:<s}, {self.empv_se[6]:>d}, {Rt_se[self.emp_se[6]]:>4.3f}, {self.i14_se[self.emp_se[6]]:>d}\n"
        i = f"#{self.emp_se[7]:<s}, {self.empv_se[7]:>d}, {Rt_se[self.emp_se[7]]:>4.3f}, {self.i14_se[self.emp_se[7]]:>d}\n"
        j = f"#{self.emp_se[8]:<s}, {self.empv_se[8]:>d}, {Rt_se[self.emp_se[8]]:>4.3f}, {self.i14_se[self.emp_se[8]]:>d}\n"
        k = f"#{self.emp_se[9]:<s}, {self.empv_se[9]:>d}, {Rt_se[self.emp_se[9]]:>4.3f}, {self.i14_se[self.emp_se[9]]:>d}\n"
        l = f"#{self.emp_se[10]:<s}, {self.empv_se[10]:>d}, {Rt_se[self.emp_se[10]]:>4.3f}, {self.i14_se[self.emp_se[10]]:>d}\n"
        m = f"#{self.emp_se[11]:<s}, {self.empv_se[11]:>d}, {Rt_se[self.emp_se[11]]:>4.3f}, {self.i14_se[self.emp_se[11]]:>d}\n"
        n = f"#{self.emp_se[12]:<s}, {self.empv_se[12]:>d}, {Rt_se[self.emp_se[12]]:>4.3f}, {self.i14_se[self.emp_se[12]]:>d}\n"
        o = f"#{self.emp_se[13]:<s}, {self.empv_se[13]:>d}, {Rt_se[self.emp_se[13]]:>4.3f}, {self.i14_se[self.emp_se[13]]:>d}\n"
        p = f"#{self.emp_se[14]:<s}, {self.empv_se[14]:>d}, {Rt_se[self.emp_se[14]]:>4.3f}, {self.i14_se[self.emp_se[14]]:>d}\n"
        q = f"#{self.emp_se[15]:<s}, {self.empv_se[15]:>d}, {Rt_se[self.emp_se[15]]:>4.3f}, {self.i14_se[self.emp_se[15]]:>d}\n"
        self.pref1 = caption + b + c + d + e
        self.pref2 = caption + f + g + h + i
        self.pref3 = caption + j + k + l + m
        self.pref4 = caption + n + o + p + q

        date_list = self.df_df['日付'].to_list()  # 日付をデータフレームからリストへ変換
        self.date_str = f"{date_list[-1]:}"  # 日付を文字列へ変換
        self.footer = '\n実効再生産数:Rt\nデータ提供NHK,独自計算,' + self.date_str
        self.mes0 = pd.Series(self.header1 + self.pref1 + self.footer)  # 地域0のメッセージ
        self.mes1 = pd.Series(self.header2 + self.pref2 + self.footer)  # 地域1のメッセージ
        self.mes2 = pd.Series(self.header3 + self.pref3 + self.footer)  # 地域2のメッセージ
        self.mes3 = pd.Series(self.header4 + self.pref4 + self.footer)  # 地域3のメッセージ
        print(self.mes0.head(), self.mes1.head(), self.mes2.head(), self.mes3.head())
        self.mes0.to_csv(self.dir + '0-------------------メッセージ4.txt', index=False, header=False,
                         sep=' ', quotechar=' ')  # 保存1
        self.mes1.to_csv(self.dir + '4-------------------メッセージ3.txt', index=False, header=False,
                         sep=' ', quotechar=' ')  # 保存2
        self.mes2.to_csv(self.dir + '8-------------------メッセージ2.txt', index=False,
                         header=False, sep=' ', quotechar=' ')  # 保存3
        self.mes3.to_csv(self.dir + '12------------------メッセージ1.txt', index=False,
                         header=False, sep=' ', quotechar=' ')  # 保存4
