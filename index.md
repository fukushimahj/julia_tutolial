# Julia チュートリアル

juliaの基本的な書き方から図の出力まで行います。

juliaの公式ドキュメントに詳しく書いてあるので、ここに記載していない事項については[こちら](https://mnru.github.io/julia-doc-ja-v1.0/index.html)で調べてください。

# Hello world!

伝統的にはじめてコードを学ぶ際には、”Hello world”と出力します。新しく”hello.jl”ファイル名で次のコードを書いて保存してください。

```julia
println("Hello world!")
```

このコードを保存した後、ターミナルでは

```bash
julia hello.jl
```

とすると実行可能です。vscodeではjuliaの拡張機能を導入している場合は、上の方にある▶️を押しても実行されます。

ここで、printlnはカッコで囲まれた結果を出力して、改行する関数となっています。lnをつけないprint()は改行なしで出力されます。次のコードも上と同じ結果を出力します。

```julia
print("Hello ")
print("world! \n")
```

# 基本的な書き方

## 変数と型

ここでは、値や文字列といった変数について説明します。juliaで扱うデータ(変数)には、データの種類を表す型が存在します。代表的なものは以下のようになります。

- 整数型(Int, Int64): 整数を表す
- 単精度浮動小数点(Float32): 4バイト実数
- 倍精度浮動小数点(Float64): 8バイト実数
- 文字列(String): 文字

juliaはpythonなどと同様に、c言語やfortranのようにコンパイル時に変数の型を宣言する必要はありませんが、型の異なる変数への代入などを行うとエラーとなる場合があります。ここで、変数の操作や計算方法を学んだ後に、異なる型への変換方法を学んでいきます。

```julia
using Printf # パッケージのimport

x = 1       # 整数
y = 1.e0    # 実数
z = 1.e3
s = "Super" # 文字列

@printf("%d %lf %le %s\n", x, y, z, s)
@printf("%s %s %s %s\n", typeof(x), typeof(y), typeof(z), typeof(s))
```

ここでは、xに整数型(Int64), yとzに実数(Float64), sに文字列(String)を代入しています。今回は、値の出力にc言語で馴染みある@printf関数を使用しています。

これを使用するために、コードの冒頭にPrintfというパッケージをusing Printfとして読み込んでいます。printfの出力方法ですが、%d, %lfは整数と実数を表し(%leは実数の指数表示), %sは文字列に対応しています。

実際に各変数の型を出力するtypeof関数を使用して、各変数の型を出力して確認しています。

また、コードにコメントを記載したい時があります。そのときは#を使用しましょう。#の後の文はコード実行時に反映されません。

## 演算

次に、変数を使用しての数値計算を行います。

```julia
using Printf
x = 5.0
y = 2.0

# 足し算, 引き算
s = x + y
v = x - y

# 掛け算, 割り算
w = x * y
r = x / y

@printf("x+y=%lf \n", s)
@printf("x-y=%lf \n", v)
@printf("x*y=%lf \n", w)
@printf("x/y=%lf \n", r)

# 整数型の割り算
a = 5
b = 2
c = a/b
println("a / b =", c)

# 自身に代入
a+=2
b*=2
@printf("a:%d b:%d \n", a, b)

# 冪
abeki = 2^2    # 4
bbeki = 10.e0^5.3
@printf("2^2=%d 10^5.3=%le \n", abeki, bbeki)

# 文字列について
s1 = "super"
s2 = "nova"
ss = s1*s2
@printf("%s \n", ss)

s1 = "Hole"
ss = "Back"*s1
@printf("%s \n", ss)

# 数学関数
sin30 = sin(pi/6.e0)
cos60 = cos(pi/3.e0)
@printf("%lf %lf\n", sin30, cos60)
```

足し算、割り算等もできます。ちなみに、整数型同士の割り算も実数型に変換してくれます(c言語等ではこうならない…)。文字列についてもつなげることができ、結合の際には*を使用します。また、数学関数であるsin()やcos()等も使用できます。

次に、型同士の変換について説明します。

**文字列から他の型**

```julia
a = parse(Int, "32")      # 整数へ
b = parse(Float64, "32")  # 実数 32.e0
c = parse(Float64, "sanjuni") # これは変換できない
```

**文字列への変換**

```julia
s12 = string(12)     # 文字列 "12" へ変換
s22 = string(22.2)   # 文字列 "22.2"へ変換
```

**実数から整数への変換**

```julia
a = Int(3.0)        # 整数3となる
b = Int(3.1)        # これはエラーとなる
c = Int(floor(3.1)) # 整数の3となる
```

# 関数

## 関数の基本的な使い型

ここでは、新しく関数を自分で作成する方法を学びます。関数は同じ処理を繰り返す時など、汎用的なコードを作成する際には必要になります。まず関数の使い方を見ていきます。

```julia
using Printf

function average(a, b)
  c = (a + b)/2
  return c
end

x = 2
y = 3
z = average(x,y)
@printf("z:%d \n", z)

```

関数はfunctionとendで囲まれた部分で構成されます。関数の名前(ここではaverage)の後に、括弧で囲まれた部分に、この関数に渡す変数(引数といいます)を書きます。今は、引数はa,bの２つしかありませんが、この数は好きな数指定できます。時には、一つも引数のない関数を定義することも可能です。関数の最後のreturn文で、結果を出力することができます。ここでは、cという変数の値を返しています。returnの後にくる変数の数も次のように好きな数指定できます。

```julia
using Printf

function tasizan_hikizan(a, b)
  c = a+b
  d = a-b
  return c, d
end

x = 2
y = 3
z, v = tasizan_hikizan(x,y)
@printf("x+y=%d, x-y=%d \n", z, v)
```

# 配列

同じ種類の変数を多数まとめて使用する場合に便利なのが配列です。また、途中でデータを加えたり、配列同士の足し算や掛け算も可能です。

## 配列の定義

配列を定義する方法はいくつかありますので、目的に応じて使い分けます。この際、変数の型を指定する必要があります。

```julia
using Printf

# 要素数0の配列を定義
x_int    = Array{Int64,1}()     # 整数型の配列
y_float  = Array{Float64, 1}()  # 実数の配列
z_moji   = Array{String, 1}()   # 文字列

# 要素数を10持つ配列を初期値0で定義
x10_int  = zeros(Int64, 10)
y10_float= zeros(Float64, 10)
```

## 要素の追加, 削除

配列の要素を追加する際には、push!()という関数を使用します。ちなみに、関数に!がつく場合は、配列の要素数を変化させるものに使用されます。削除する方法もあり、ここではdeleteat!()を使用する方法を紹介します。

```julia
# 要素数0の配列を定義
x_int    = Array{Int64,1}()     # 整数型の配列
y_float  = Array{Float64, 1}()  # 実数の配列
z_moji   = Array{String, 1}()   # 文字列

# 要素を追加
push!(x_int, 1)
push!(x_int, 2)
push!(x_int, 3)
println(x_int)

push!(y_float, 868.e0)
push!(y_float, 657.e0)
push!(y_float, 567.e0)
println(y_float)

push!(z_moji, "Ikeda")
push!(z_moji, "Fujioka")
push!(z_moji, "Nakamura")
println(z_moji)

# 要素を削除
deleteat!(x_int, 1)             # 配列の1番目を削除(1)
println("x_int:", x_int)

deleteat!(y_float, 2)           # 配列の２番目を削除(657.e0)
println("y_float:", y_float)

deleteat!(z_moji, 3)            # 配列の3番目を削除("Nakamura")
println("z_moji:", z_moji)

# 配列の要素数を取得
println(size(x_int)[1])
```

ちなみに、配列のindexについてですが、c言語やpythonではindexは０から始まりますが、juliaでは1から始まるので注意が必要です(fortranもそうです)。

## 配列の結合

2つの配列を結合する場合は次のようにします
```julia
a = [1,2,3]
b = [4,5,6,7]

append!(a, b) # aにbの要素をつけ加える

println(a)
```


## 多次元配列

多次元配列の定義も同様にできます。２列目以降を追加する際には;で区切る。

```julia

matrixA = [ 1 2; 3 4]   # 2 X 2の行列
println(matrixA)

matrix0 = zeros(Float64, 2, 2) # 最初に0を持つ0x0の行列
```

## 配列同士の演算

配列同士の各要素ごとについて計算することも可能です。この場合、演算子の前に.を追加する必要があります。

```julia
# 配列を定義
a = [ 1, 2, 3, 4]
b = [ 3, 6, 9, 12]

# 行列の各要素ごとに足し算
sum_ab = a .+ b 
println(sum_ab)

# 掛け算
time_ab = a .* b
println(time_ab)

# 定数をかける
a2 = 2 .* a
println(a2)

# 冪乗
aa = a .^2
println(aa)

# 文字列への変換
string_a = string.(a)
println(string_a)

# 文字列から整数へ
int_a = parse.(Int64, string_a)
println(int_a)
```
配列は型が同じものしか格納できない制限があります。異なる型の配列を含みたい場合などに次の辞書式が有効に使えます。

## 辞書式

配列は要素番号に対して値が振り分けられていますが、辞書式は任意の数字、もしくは文字を配列として格納することができます。辞書式はDict{,}()で定義します。最初の文字で、indexとして使用する数値、もしくは文字、次の部分にそれが示す変数の型もしくは配列(辞書式)を指定できます。実際に見てみましょう。

```julia
# 辞書式

m = Dict{String, Int64}() # 文字式をindexに指定し、整数型を格納
m["Ojima"]    = 14
m["Taneichi"] = 16
m["Sasaki"]   = 17
println(m)     # 文字式全体を出力
println("Taneichi:", m["Taneichi"]) # 各indexに格納した値を出力

# m = Dict{Int64, Float64}などindexには整数等も使用かのうです。

# 辞書式の要素に配列を指定する場合
a1 = [1,2,3,4,5]
a2 = [2,4,6,8,10]
dd = Dict{Int64, Array{Int64, 1}}()    # 整数型の配列を持つ辞書を定義
dd[1] = a1
dd[2] = a2
println(dd[1], dd[2])                  
println(dd[1][2])        # 各要素を出力することも可能です。

# 異なる型を持つ辞書も定義することも可能です。
d2 = Dict{Int64, Any}()  #この場合はanyを指定します。
d2[1] = m      # 要素として辞書を指定することも可能です。
d2[2] = a1
d2[3] = 3.14
println(d2)

# 辞書を要素に持つ辞書も定義可能です
d3 = Dict{Int64, Dict{Int64, Float64}}()

```

# 条件式

変数がある条件を満たす場合だけ計算がしたい場合など、if文が非常に有効となります。if文はそれに続く条件式が満たされる(true)場合だけif文内の処理が実行されます。また、満たされない(false)場合についてもelse文を用いることで、その場合の処理も指定できます。

```julia
x = 10

if x > 1
  println("x > 1") # 条件が満たされるのでこの文が実行される 
end

if x > 10
  println("x > 10") # 条件が満たされないので、ここの部分は実行されない
else
  println("x <= 10") # if文の条件がfalseなので、こちらが実行される。
end

y = 20
if (x > 1 && y > 10)   # 条件式を2つ指定することも可能で、&& はどちらもtrueの場合だけtrueとなります。
  println("x > 1 && y > 10 is true")
end

if (x > 1 || y < 10)   # どちらか(もしくは両方)がtrueの場合に実行したい場合は||を使用します　
　println("x > 1 || y < 10 is true")
end
```

# 繰り返し処理

プログラムでは、同じ処理を何度も繰り返した場合があります。この場合には、for文やwhile文などを使用します。

```julia
# for文

sum = 0
a   = Array{Int64, 1}() 
a2  = zeros(Int64, 10)
for i=1:10
  global sum += i # globalが必要ですが、これは後のスコープで説明します。 
  push!(a, i)     # 各iを配列に加えていきます。
  a2[i] = 2*i     # あらかじめ用意した配列の要素にも代入できます。
end

println("sum:", sum)
println("a:", a)
println("a2:", a2)

# 配列の各要素についてのループ
for aa in a2      # ここでは、a2に格納された10この数字が順番にaaへと代入されます。
  println(aa)
end

# 文字列の配列についても可能です
mr = ["sasaki", "ojima", "taneichi", "nishino", "mercedes", "mori"]
for sp in mr
  println(sp)
end

```

for文では、i=1:10となっている場合、iに順番に1,2,…,10を代入した場合の計算が行われます。例えば、1から10までの足し算を行うなど一つずつ記載するのは大変ですが、このように簡単化できます。また、配列の各要素についてもループすることができます。

for文を途中で抜け出したい場合などは、break文を使用します。

```julia
for i=1:100
  println(i)
  if(i == 10) # if文は次の章で説明します。...
    break     # ここでループを抜ける
  end
end
```

次にwhile文について説明します。while文では、条件式が真の間は制限なく繰り返し処理が行われます。

```julia

i   = 0
while i<10 # i<10が条件式なので、iが9まで処理が繰り返されます。
  println(i)
  global i+=1
end

i = 100
while i<10 # この場合は、whileの最初から条件式が否なので、while文の中の処理は一度も実行されません。
  println("!!!!!!!")
end
```

while文は, 時間発展を計算する際に、ある時刻まで方程式を積分したい場合などで非常に有効です。

# スコープ

juliaの機能の一つに変数のスコープ(有効範囲)があります。変数には型の他に、スコープという属性がつけられます。スコープにはlocalとglobalがあり、複雑な振る舞いをする場合があります。例えば、次のコードを見てください。

```julia
x = 1  
for i=1:10
  x = 10
end
println(x)
```

このコードは一見するとループの中で、xに10を代入しているので、結果も10が出力されますが、結果としては１が出力されます。これは、冒頭でx=1と定義されたxはglobal変数であるのに対し、for文の中でx=10で定義されているものはfor文の中でのみ有効なlocal変数となります。この二つの変数は、xで定義されていますが、別の変数として扱われます。

では、for文の中でglobal変数のxについて計算したい場合は、どうすればいいでしょうか。その場合には、globalを変数の前につけます。

```julia
x = 1   # global 変数
for i=1:10
  x = 10 # これはlocal変数のx
  global x = 20 # 変数の前にglobalをつけるとglobal変数について計算できます
end
println(x)
```

結果として、20が出力されて、for文の中の計算が反映されています。

ちなみに、local変数はfor やfunctionのendまでしか有効ではなく、それ以降は使用できません。例えば、次のようなコードはエラーとなります。

```julia
for i=1:10
  y = 30    # local 変数
end
println(y)  # これは使用できない...
```

## スコープを避ける…

上の説明を見て？となったり、いちいち考えるのが面倒だと思った人も多いと思います(わたしもスコープは苦手です)。なので、計算の主な部分は関数の中で行うことをお勧めします。関数の中で宣言した変数は全てlocal変数となるので、上記のような問題は起きなくなります。

```julia
function test1()
  x = 1             # function test1の local 変数として定義  
  for i=1:10
    x = 10　　　     # for 文内部で新たにlocal変数として認識されることなく、これもfunction test1のlocal変数
    global y = 20   # global属性がついていないと、function test1では宣言されていないので、for文のlocal変数として定義されてしまう。このため、エラーが出力される
  end
  println(x)
  println(y)
end
test1()
println(y)          # global属性なので、関数の外側でも使用可能
```

最初に見せたように、関数内部のlocal変数としてxは宣言しているので、for文の中では新たにlocal変数としては定義されず、そのまま値を代入することができます。一方、for文の中であらたに値を定義した場合は、その外では使用できなくなります。

ややこしい概念ですが、基本for文を使用するときは、関数を定義して使うくらいでもいいかもしれません。

ちなみに、配列はそのままglobal変数のindexを指定しているだけなので、for文内でもglobal変数にそのまま代入できます。ただし、for文内部で新たに定義すると、local変数となります。

```julia
x = zeros(Int64, 10)
for i=1:10
  x[i] = i   # これはglobal変数
end
println(x)

for i=1:10
  x 　　= zeros(Int64, 10)  # ここでlocal変数を宣言
  x[i] = 33                # ここもlocal変数
end
println(x)                   # 上のループの変更は反映されない
```

# パッケージ取得

julia では様々パッケージを使用することができ、図を描くなども簡単にできようになっています。ターミナルを使用している人は、ターミナル上でjuliaとコマンドし、julia>の次のように入力してみてください。(VSCODEの人は、再生ボタンを押すと同じものがでてくると思います。)

```julia
import Pkg; Pkg.add("PyPlot")
```

ここでは、”PyPlot”という図を描く際に使用するパッケージをインストールしました。

# ファイル読み書き

ここではファイルへの書き出し、読み込み方法を紹介します。ここでの方法は一例となります。

```julia
using Printf

ns = 0:0.01:2
sin_s = sin.(ns.*pi)   # [0, 2pi]でsin関数を準備

f = open("./data.dat", "w")

# ファイルへの書き込み
for (i, n) in enumerate(ns) # enumerateは配列の要素順に、i=1,2,3..., n=ns[1], ns[2], ns[3], ...と代入してくれる関数です。
  @printf(f, "%d %e %e \n", i, n, sin_s[i])
end

close(f) # closeでファイルの書き込み終わり
```

openという関数では、データを書き込みし(読み込み)たいファイルの名前を指定します。”w” (”r”) はデータを書き込み(読み込み)することを表しています。

ここで、作成した./data.datというファイルを今度は読み込みしてみましょう。

```julia
using Printf

# 各行ごとに読み込み(linesに文字列配列として格納)　
lines = open("./data.dat", "r") do fp
  readlines(fp)
end

for line in lines
  d = parse.(Float64, split(line)) # splitで空白ごとに、配列として分割 例えば, "ibaraki tsukuba" => ["ibaraki", "tsukuba"], parse.で一括で実数に変換
  i     = Int(d[1])  # 配列のindex
  theta = d[2]*pi    # theta
  sin   = d[3]	     # sin theta
  @printf("%d %e %e \n", i, theta, sin)
end                                 
																	 
```

# 図の作成

図はPyPlotというパッケージを使用して作成してみましょう。ここでは、試しに上で作成したファイルを使用して、図を作成してみます。

```julia
using PyPlot

# 各行ごとに読み込み(linesに文字列配列として格納)　
lines = open("./data.dat", "r") do fp
  readlines(fp)
end

thetas = Array{Float64, 1}()
sins   = Array{Float64, 1}()
for line in lines
  d = parse.(Float64, split(line)) # splitで空白ごとに、配列として分割 例えば, "ibaraki tsukuba" => ["ibaraki", "tsukuba"], parse.で一括で実数に変換
  i     = Int(d[1])  # 配列のindex
  theta = d[2]*pi    # theta
  sin   = d[3]	     # sin theta

  push!(thetas, theta)
  push!(sins, sin)
end                                 
						
figure()
plot(thetas, sins)		
show()									 
```

上では、sin関数を描いています。

図はlabelや色の指定など様々なものが設定できるので、試してみましょう。また、図を表示するだけではなく、保存することも可能です。

```julia
using PyPlot

# 各行ごとに読み込み(linesに文字列配列として格納)　
lines = open("./data.dat", "r") do fp
	readlines(fp)
end

thetas = Array{Float64, 1}()
sins   = Array{Float64, 1}()
coss   = Array{Float64, 1}()
for line in lines
  d = parse.(Float64, split(line)) # splitで空白ごとに、配列として分割 例えば, "ibaraki tsukuba" => ["ibaraki", "tsukuba"], parse.で一括で実数に変換
  i     = Int(d[1])  # 配列のindex
  theta = d[2]*pi    # theta
  sin   = d[3]	     # sin theta

  push!(thetas, theta)
  push!(sins, sin)
  push!(coss, cos(theta))
end                                 
						
# 図の作成

figure(figsize=(8,6))     # figsizeでサイズも指定可能

#plot関数で実際に図を描く
plot(thetas, sins, "-" , lw=2, label=L"\sin") 		
plot(thetas, coss, "--", lw=2, label=L"\cos") 	

# "-", "--"は実線か破線の線の種類の選択
# lw: 線の太さ
# labelは次のlegend関数が使用された際のラベルの設定
# 他にもcolorなどもある

legend(loc="lower left", fontsize=18) # fontsize は文字の大きさ
xlabel(L"\theta", fontsize=18)        # xlabel, ylabelを設定
ylabel(L"\sin(\theta) \, \cos(\theta)", fontsize=18)

xlim(0, 2*pi) # x軸 & y軸の範囲を指定
ylim(-1,1)
					 
savefig("./sincos.png") # 図の名前が必要、png, jpeg, pdfと好きな形式を選べる
show()				
close()
```
