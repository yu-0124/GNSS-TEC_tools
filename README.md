# rdrnxの使用法

RINEX形式のデータファイル(AAAADDD0.YYo AAAAは局ID、DDDは通算日, YYは年の下二桁)を読み込んでL4 (=L1-L2)を作成して衛星毎にソートして出力する。単位は時刻がhour (UT)、L4はTECU (TEC unit, 10E16 electrons per square meter)である。

使用例：

　　rdrnx < 04932030.09o 

RINEX fileである04932030.09o (屋久島の950493局の2009年7月22日、通算日で203日)のGPSデータファイルを読んでL4を衛星毎に出力する。

出力例： 各行は時刻（UT, hour）とL4 (TECU)から成る。衛星毎に>で始まるヘッダーが付いており、衛星番号、データ数、局名、除去した下駄履きの値、が示されている。

```
> sat# 9 #data: 388 site:0493 bias:    -48954.0352
 20.767     -0.0006
 20.775      0.1273
 20.783      0.1253
 20.792      0.1556
 20.800      0.0724
 20.808      0.2072
 20.817      0.3060
 20.825      0.5180
```

# rdephの使用法

RINEX形式の軌道ファイル(AAAADDD0.YYn AAAAは局ID、DDDは通算日、YYは年の下二桁）を読み込んで、地球固定座標系における各衛星の位置（XYZ座標値、単位m）を３分(0.05時間）毎に出力する。

使用例：

　　rdeph < 04932030.09n

RINEX軌道ファイルである04932030.09n（屋久島の950493局の2009年7月22日、通算日で203日）のGPS軌道ファイルを読み込んで衛星の位置のXYZ値を３分毎に出力する。

出力例：
```
    0.00  2 -.13716196699E+08 -.10694036001E+08 -.20257051753E+08
    0.00  3 0.24392166998E+08 0.83771893845E+07 0.67596149148E+07
    0.00  4 -.77670740916E+07 -.21002870036E+08 -.14133298063E+08
    0.00  6 0.23138749340E+08 0.12740274927E+08 0.33411813835E+07
    0.00  7 0.70645039310E+07 -.24861868739E+08 -.58263967196E+07
    0.00  8 -.26397518768E+06 -.25041812066E+08 0.79753308477E+07
    0.00  9 -.13805186266E+08 0.12232787608E+08 0.18418895377E+08
    0.00 10 -.23049147807E+08 -.11157741823E+07 -.13603407392E+08
    0.00 11 0.13101182698E+08 -.14350313030E+08 0.17756192825E+08
    0.00 12 -.23487165715E+08 0.11921537748E+08 -.27720430004E+07
```

最初の数が時刻(UT, hour)、二番目の整数が衛星番号、その後に続く三つの実数が衛星の地球固定系における位置のXYZ座標値(meter)である。

# rdrnx3の使用法

バージョン３のRINEXファイル名、衛星システムを表すアルファベット１文字、さらに取り出す時間帯の最初と最後の時間（ＵＴ、単位は時間）を引数で与える。アルファベットはGPSがG,GLONASSがR, GalileoがE, QZSSがJといった具合である。

使用例：

　　rdrnx3 04932030.20o J 0.0 24.0 

Multi-GNSSのRINEXファイル(v.3)を読んでQZSSのデータを一日分(0時から24時まで)L4に変換して衛星ごとにソートして出力する。

# tomo （三次元トモグラフィー） の使用法

下記論文を参照のこと
- He, L. and K. Heki, Three-dimensional tomography of ionospheric anomalies immediately before the 2015 Illapel earthquake, central Chile, J. Geophys. Res. Space Phys., 123, 4015-4025, doi:10.1029/2017JA024871, 2018.
- Muafiry, I. N., K. Heki, and J. Maeda, 3D tomography of mid-latitude sporadic-E in Japan from GNSS-TEC data, Earth Planets Space,70, 45, doi:10.1186/s40623-018-0815-7, 2018.
- He, L., K. Heki,and L. Wu, Three-dimensional and trans-hemispheric ionospheric electron density changes by the great Solar eclipse in North America on 21 August 2017, Geophys. Res. Lett., 45, 10,933-10,940, doi:10.1029/2018GL080365, 2018.
- Muafiry, I.N.and K. Heki, 3D tomography of the ionospheric anomalies immediately before and after the 2011 Tohoku-oki (Mw9.0) earthquake, J. Geophys. Res. Space Phys., in press.
