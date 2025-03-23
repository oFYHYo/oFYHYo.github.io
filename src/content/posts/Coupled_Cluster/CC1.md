---
title: Coupled Cluster理论学习笔记（1）：耦合簇形式理论
published: 2025-03-23
description: ''
image: ''
tags: [电子结构]
category: '电子结构'
draft: false 
lang: 'zh_CN'
---

&emsp;&emsp;最近在写CC。因为CC公式大多比较长，而且推起来很麻烦，所以我就打算整理一下公式再写几篇博客，以防之后忘记。

&emsp;&emsp;Coupled Cluster（耦合簇，CC）理论是一种量子化学计算方法，它能够恰当地考虑电子的动态相关，从而能够高精度地计算一些体系的相关能。CC 最早被核物理学家所提出，主要用来计算核反应体系的一些性质[^1]。上世纪六、七十年代，化学家们将CC用于处理原子分子体系，开发了CCSD、CCSD（T）等一系列计算方法，并且取得了巨大的成就。时至至今日，CCSD（T）仍被认为是单参考态方法的金标准，被广泛应用于小分子体系的精确计算。

&emsp;&emsp;作为本系列博客的第一篇文章，本文将会介绍CC形式理论，并会给出CC能量与振幅所满足的一般方程。在阅读本系列文章之前，建议读者至少要有一定的量子力学基础。

## 1. 组态波函数的二次量子化表示

### 1.1 精确波函数的 CI（Configure Interaction）展开

&emsp;&emsp;对于一特定的 n 电子波函数，可用符号 $\ket{p,q,\cdots r,s}$ 表示，该符号代表第一个电子占据第 $p$ 个单电子态（一般取为HF分子轨道或自旋轨道）、第二个电子占据第 $q$ 个单电子态，以此类推。像这样电子在单电子态上的排列被称为组态。所有可能的组态波函数与真空态所构成的态矢量集可以张成一个空间，该空间一般被称为Fock空间（或组态空间）。在Fock空间中，所有态矢量都能用组态波函数的线性组合来表示，即

$$
\begin{equation}
\ket{\Psi} = \sum C_{pq\cdots rs}\ket{p,q,\cdots r,s}
\end{equation}
$$

如果我们选定一个态作为参考态（如 HF 基态），把其他态看作参考态的激发，上式就可以写成（在本文中，$ijkl$ 指代占据轨道，$abcd$ 指代非占据轨道，$pqrs$ 指代任意轨道）

$$
\begin{equation}
\ket{\Psi} = \ket{\Psi_0}+\sum_{ia}c_i^a\ket{\Psi_i^a}+\sum_{iajb}c_{ij}^{ab}\ket{\Psi_{ij}^{ab}}+\cdots
\end{equation}
$$

其中 $\ket{\Psi_0}$ 代表参考态，$\ket{\Psi_i^a}$ 代表将参考态位于第 $i$ 个单电子态上的电子给激发到第 $a$ 个单电子态，以此类推，直至遍历所有可能的组态。这便是精确波函数的 **CI 展开**。

### 1.2 组态波函数的二次量子化表示

&emsp;&emsp;写出组态波函数二次量子化表示

$$
\begin{equation}
\ket{p,q,\cdots r,s} = \hat{a}^\dagger_s\hat{a}^\dagger_r\cdots\hat{a}^\dagger_q\hat{a}^\dagger_p\ket{}
\end{equation}
$$

式中 $\hat{a}^\dagger$ 和 $\hat{a}$ 分别代表费米子的产生算符和湮灭算符，它们之间所满足的对易关系如下

$$
\begin{equation}
\begin{aligned}
[\hat{a}^\dagger_p,\hat{a}^\dagger_q]_+ &= \hat{a}^\dagger_p\hat{a}^\dagger_q+\hat{a}^\dagger_q\hat{a}^\dagger_p=0\\
[\hat{a}_p,\hat{a}_q]_+ &= \hat{a}_p\hat{a}+\hat{a}_q\hat{a}_p=0\\
[\hat{a}^\dagger_p,\hat{a}_q]_+ &= \hat{a}^\dagger_p\hat{a}_q+\hat{a}_q\hat{a}^\dagger_p=\delta_{pq}
\end{aligned}
\end{equation}
$$

通过上述对易关系，我们可以验证一下费米子体系的泡利不相容原理

$$
\begin{equation}
\ket{p,q,\cdots r,s} = \hat{a}^\dagger_s\hat{a}^\dagger_r\cdots\hat{a}^\dagger_q\hat{a}^\dagger_p\ket{} = -\hat{a}^\dagger_s\hat{a}^\dagger_r\cdots\hat{a}^\dagger_p\hat{a}^\dagger_q\ket{} = \ket{q,p,\cdots r,s}
\end{equation}
$$

即交换两电子位置，波函数变号。也可以通过产生湮灭算符的对易关系和真空态的一些性质，来计算一些积分，比如用于计算两个双电子态的重叠积分

$$
\begin{equation}
\begin{aligned}
\braket{p,q|r,s} &= \bra{}\hat{a}_p\hat{a}_q\hat{a}^\dagger_s\hat{a}^\dagger_r\ket{} = \bra{}\hat{a}_p(\delta_{qs}-\hat{a}^\dagger_s\hat{a}_q)\hat{a}^\dagger_r\ket{}\\
&= \delta_{qs}\bra{}\hat{a}_p\hat{a}^\dagger_r\ket{}-\bra{}\hat{a}_p\hat{a}^\dagger_s\hat{a}_q\hat{a}^\dagger_r\ket{}\\
&= \delta_{pr}\delta_{qs}-\delta_{qr}\bra{}\hat{a}_p\hat{a}^\dagger_s\ket{}\\
&= \delta_{pr}\delta_{qs}-\delta_{ps}\delta_{qr}
\end{aligned}
\end{equation}
$$

这就是双电子情况下不同组态之间的正交归一关系。通过二次量子化语言，我们可以在不写出组态波函数的显式表达式的情况下，仅用代数方法便可求得各种积分。但对于实际体系，由于电子数非常多，直接使用（6）式中所示方法就需要处理大量产生湮灭算符的乘积，这将会变得非常困难，有没有方法能够避免这个问题呢？

## 2. 费米真空与正则序

### 2.1 费米真空与 q-产生湮灭算符

&emsp;&emsp;答案当然是有的。回顾波函数的 CI 展开式（2），我们可以把任意组态看作参考态的激发，即将参考态看作真空态，将其他组态写成产生湮灭算符作用在参考态上的形式，即

$$
\begin{equation}
\begin{aligned}
\ket{\Psi_0} &= \ket{\Psi_0} \\
\ket{\Psi_i^a} &= \hat{a}^\dagger_a\hat{a}_i\ket{\Psi_0}\\
\ket{\Psi_{ij}^{ab}} &= \hat{a}^\dagger_a\hat{a}^\dagger_b\hat{a}_i\hat{a}_j\ket{\Psi_0}\\
& \ldots
\end{aligned}
\end{equation}
$$

通过把真空态与一部分产生算符合并成一个参考态，我们可以极大地减少要处理的产生湮灭算符的个数。如果选取的参考态为 HF 基态，我们就称该参考态为**费米真空**。

&emsp;&emsp;在将真空态换为费米真空后，产生湮灭算符的一些性质也会发生些许改变。原本湮灭算符作用在真空态上会等于零，即  $\hat{a}_p\ket{}=0$，但是在把真空态换成费米真空后，这个等式并不一定成立：非占据轨道的湮灭算符作用在 HF 基态上等于零，而占据轨道的湮灭算符作用在 HF 基态上并不为零，相反，占据轨道的产生算符作用到 HF 基态后会变成零，即

$$
\begin{equation}
\begin{aligned}
\hat{a}_a\ket{\Psi_0}=0\\
\hat{a}^\dagger_i\ket{\Psi_0}=0
\end{aligned}
\end{equation}
$$

可以看到，在选择费米真空作为真空态后，占据轨道的产生算符会表现出‘湮灭算符‘的性质，占据轨道的湮灭算符会表现出‘产生算符’的性质。因此，在选择费米真空作为真空态后，我们可以将占据轨道的产生算符与非占据轨道的湮灭算符归为一类‘湮灭算符’、将占据轨道的湮灭算符与非占据轨道的产生算符归为一类‘产生算符’。

&emsp;&emsp;为了与通常意义下的产生湮灭算符做区别，我们通常会把占据轨道的产生算符与非占据轨道的湮灭算符称为 **$q$-湮灭算符**、将占据轨道的湮灭算符与非占据轨道的产生算符称为 **$q$-产生算符**。与产生湮灭算符的对易关系相对应，我们也可以求出 $q$-产生湮灭算符之间的对易关系，即

$$
\begin{equation}
\begin{aligned}
[\hat{a}^\dagger_p,\hat{a}^\dagger_q]_+ &= \hat{a}^\dagger_p\hat{a}^\dagger_q+\hat{a}^\dagger_q\hat{a}^\dagger_p=0\\
[\hat{a}_p,\hat{a}_q]_+ &= \hat{a}_p\hat{a}+\hat{a}_q\hat{a}_p=0\\
[\hat{a}^\dagger_i,\hat{a}_j]_+ &= \hat{a}^\dagger_i\hat{a}_j+\hat{a}_j\hat{a}^\dagger_i=\delta_{ij}\\
[\hat{a}^\dagger_a,\hat{a}_b]_+ &= \hat{a}^\dagger_a\hat{a}_b+\hat{a}_b\hat{a}^\dagger_a=\delta_{ab}\\
[\hat{a}^\dagger_i,\hat{a}_a]_+ &= \hat{a}^\dagger_i\hat{a}_a+\hat{a}_a\hat{a}^\dagger_i=0\\
[\hat{a}^\dagger_a,\hat{a}_i]_+ &= \hat{a}^\dagger_a\hat{a}_i+\hat{a}_i\hat{a}^\dagger_a=0
\end{aligned}
\end{equation}
$$

&emsp;&emsp;下面我们来讨论 $q$-产生湮灭算符的物理意义。我们知道，通常的产生（湮灭）算符作用到一个态上后会在指定电子态上增加（减少）一个电子。在定义费米真空后，我们将原本的电子态划分成了占据与非占据两种电子态，对于占据轨道的产生湮灭算符，其物理意义没变：将其作用到一个态上后会在非占据轨道上增加（减少）一个电子；但对于占据轨道的产生湮灭算符，可以认为在占据轨道中湮灭一个电子就相当于产生一个“空穴”，在空穴中产生一个电子就相当于湮灭了一个“空穴”，所以占据轨道的产生（湮灭）算符会具有 $q$-湮灭（产生）算符的性质，占据轨道的产生湮灭算符实际上对应着空穴的产生湮灭而非电子的产生湮灭。综上所述，**$q$-产生算符可以产生电子或空穴、$q$-湮灭算符可以湮灭电子或空穴**，这种引入电子和空穴的概念的表述也被称为 **Partical-Hole 表述**。

&emsp;&emsp;以一张表来归纳真空与费米真空之间的异同

||真空|费米真空|
|:----:|:----:|:----:|
|真空态|$\ket{}$|$\ket{\Psi_0}$|
|产生算符|$\hat{a}^\dagger_p$|$\hat{a}^\dagger_a,\hat{a}_i$|
|湮灭算符|$\hat{a}_p$|$\hat{a}^\dagger_i,\hat{a}_a$|
|产湮算符的物理意义|产生（湮灭）电子|占据轨道的算符为产生（湮灭）空穴，非占据轨道为产生（湮灭）电子|
|非零对易关系|$[\hat{a}^\dagger_p,\hat{a}_q]_+ =\delta_{pq}$|$\begin{aligned}[\hat{a}^\dagger_i,\hat{a}_j]_+ = \delta_{ij},[\hat{a}^\dagger_a,\hat{a}_b]_+ = \delta_{ab}\end{aligned}$|

&emsp;&emsp;用一个例题来加深对这些概念的理解

**例1.** 令 $\hat{H} = \sum h_{pq}\hat{a}^\dagger_p\hat{a}_q$，计算算符 $\hat{H}$ 在 HF 基态下的期望。

答:

$$
\begin{equation}
\begin{aligned}
\bra{\Psi_0}\hat{H}\ket{\Psi_0} &= \bra{\Psi_0}\sum h_{pq}\hat{a}^\dagger_p\hat{a}_q\ket{\Psi_0}\\
&=\sum h_{pq} \bra{\Psi_0}\hat{a}^\dagger_p\hat{a}_q\ket{\Psi_0}\\
&=\sum h_{ij} \bra{\Psi_0}\hat{a}^\dagger_i\hat{a}_j\ket{\Psi_0}+\sum h_{ia} \bra{\Psi_0}\hat{a}^\dagger_i\hat{a}_a\ket{\Psi_0}\\
&+\sum h_{ai} \bra{\Psi_0}\hat{a}^\dagger_a\hat{a}_i\ket{\Psi_0}+\sum h_{ab} \bra{\Psi_0}\hat{a}^\dagger_a\hat{a}_b\ket{\Psi_0}
\end{aligned}
\end{equation}
$$

其中

$$
\begin{equation}
\begin{aligned}
&\bra{\Psi_0}\hat{a}^\dagger_i\hat{a}_j\ket{\Psi_0}=\bra{\Psi_0}\delta_{ij}-\hat{a}_j\hat{a}^\dagger_i\ket{\Psi_0}=\delta_{ij}\\
&\bra{\Psi_0}\hat{a}^\dagger_i\hat{a}_a\ket{\Psi_0}=0\\
&\bra{\Psi_0}\hat{a}^\dagger_a\hat{a}_i\ket{\Psi_0}=-\bra{\Psi_0}\hat{a}_i\hat{a}^\dagger_a\ket{\Psi_0}=0\\
&\bra{\Psi_0}\hat{a}^\dagger_a\hat{a}_b\ket{\Psi_0}=0
\end{aligned}
\end{equation}
$$

将（11）代入（10）中，可得

$$
\begin{equation}
\bra{\Psi_0}\hat{H}\ket{\Psi_0} = \sum h_{ij}\delta_{ij}=\sum h_{ii}
\end{equation}
$$

### 2.2 算符的正则序与wick定理

&emsp;&emsp;从（6）式和（11）式的计算过程来看，无论是真空还是费米真空，计算这类含产生湮灭算符的积分的主要步骤实际上是通过对易关系，将原本的产生湮灭算符的乘积给化成左边是产生算符右边是湮灭算符的形式，容易证明，这种左产生右湮灭形式的算符作用到真空态上会等于零。我们将这种左产生又湮灭形式的称为算符的 **正则序（normal-ordering）** 形式，定义正则序算符，其作用是将作用到的算符变换为左产生右湮灭的形式，用符号 $\{\cdots\}_v$ 来表示，即

$$
\begin{equation}
\{\hat{X}\hat{Y}\hat{Z}\cdots\hat{U}\hat{V}\hat{W}\}_v = (-1)^P(\underbrace{\hat{U}\hat{X}\cdots\hat{Y}}_{\text{产生算符}}\underbrace{\hat{V}\cdots\hat{W}\hat{Z}}_{\text{湮灭算符}})
\end{equation}
$$

上式中，$P$ 等于从原算符乘积的排列到算符的正则序的排列对换的次数，即置换的逆序数。容易证明，算符的正则序作用到真空态后会等于零。通过引入正则序，我们就可以将形如例1的算符期望的计算问题给转化成通过对易关系，将原算符给化成正则序形式的算符，而wick定理就能做到这一点。

&emsp;&emsp;在介绍wick定理之前，我们先介绍一个概念——**算符的收缩**。两个算符 $\hat{A},\hat{B}$ 的收缩可以由符号 $\hat{A}^\bullet\hat{B}^\bullet$ 来表示[^2]。算符收缩的定义式为

$$
\begin{equation}
\hat{A}^\bullet\hat{B}^\bullet=\hat{A}\hat{B}-\{\hat{A}\hat{B}\}_v
\end{equation}
$$

即一个算符的收缩等于算符本身减去算符的正则序。根据正则序的定义式（13）和收缩的定义式（14），我们可以计算后面需要用到的 $q$-产生湮灭算符之间的收缩。$q$-产生湮灭算符的收缩的非零项如下

$$
\begin{equation}
\begin{aligned}
\hat{a}^{\dagger\bullet}_i\hat{a}_j^\bullet = \delta_{ij}\\
\hat{a}^{\bullet}_a\hat{a}_b^{\dagger\bullet} = \delta_{ab}
\end{aligned}
\end{equation}
$$

在定义了算符的收缩后，我们就能够介绍了。wick定理可以表述为：一个算符可以等于该算符的正则序与对该算符所有可能的收缩的正则序之和，即

$$
\begin{equation}
\begin{aligned}
\hat{A}\hat{B}\cdots\hat{C}\hat{D} &= \{\hat{A}\hat{B}\cdots\hat{C}\hat{D}\}_v\\
&+\{\hat{A}^\bullet\hat{B}\cdots\hat{C}^\bullet\hat{D}\}_v+\text{所有可能的单收缩的正则序}\\
&+\{\hat{A}^\bullet\hat{B}^{\bullet\bullet}\cdots\hat{C}^\bullet\hat{D}^{\bullet\bullet}\}_v+\text{所有可能的双收缩的正则序}\\
&+\cdots
\end{aligned}
\end{equation}
$$

wick定理的一个简单证明过程可参见[^3]。通过wick定理，我们可以很方便地将原本的多算符乘积化成一系列正则序的求和，从而简化一些计算。在之和推导CC能量和振幅的具体表达式的时候会经常用到这个定理。

## 3. CC形式理论

### 3.1 精确波函数的簇展开

&emsp;&emsp;下面我们将从波函数的 CI 展开出发，逐步引入耦合簇理论的一些基本概念。

&emsp;&emsp;使用上一节引入的费米真空与产生湮灭算符，我们可以把精确波函数的 CI 展开式（2）写成一个激发算符作用在 HF 基态上的形式，即

$$
\begin{equation}
\begin{aligned}
\ket{\Psi} &= \ket{\Psi_0}+\sum_{ia} c_i^a\hat{a}^\dagger_a\hat{a}_i\ket{\Psi_0}+\sum_{iajb}c_{ij}^{ab}\hat{a}^\dagger_a\hat{a}^\dagger_b\hat{a}_i\hat{a}_j\ket{\Psi_0}+\cdots\\
&= (1+\sum_{ia} c_i^a\hat{a}^\dagger_a\hat{a}_i+\sum_{iajb}c_{ij}^{ab}\hat{a}^\dagger_a\hat{a}^\dagger_b\hat{a}_i\hat{a}_j+\cdots)\ket{\Psi_0}\\
&= \hat{C}\ket{\Psi_0}
\end{aligned}
\end{equation}
$$

其中

$$
\begin{equation}
\begin{aligned}
\hat{C}=1+\sum_{ia} c_i^a\hat{a}^\dagger_a\hat{a}_i+\sum_{iajb}c_{ij}^{ab}\hat{a}^\dagger_a\hat{a}^\dagger_b\hat{a}_i\hat{a}_j+\cdots
\end{aligned}
\end{equation}
$$

这就是 CI 展开的二次量子化表示。引入簇算符 $\hat{T}$，其定义为

$$
\begin{equation}
\hat{T} = \sum_{n=1}^\infty\hat{T}_n
\end{equation}
$$

其中

$$
\begin{equation}
\begin{aligned}
\hat{T}_1 &= \sum t_i^a\hat{a}^\dagger_a\hat{a}_i\\
\hat{T}_2 &= \sum t_{ij}^{ab}\hat{a}^\dagger_a\hat{a}^\dagger_b\hat{a}_i\hat{a}_j\\
&\cdots
\end{aligned}
\end{equation}
$$

其中各项系数一般被称为 CC **振幅**。令

$$
\begin{equation}
\hat{C} = 1+\hat{T}+\frac{\hat{T}^2}{2!}+\cdots=e^{\hat{T}}
\end{equation}
$$

可以由上式计算出算符 $\hat{C}$ 中系数与算符 $\hat{T}$ 中系数之间的对应关系，将（21）式与（18）式相比较，可得

$$
\begin{equation}
\begin{aligned}
c_i^a = &t_i^a\\
c_{ij}^{ab} = &t_{ij}^{ab}+t_i^at_j^b-t_i^bt_j^a\\
c_{ijk}^{abc} = &t_{ijk}^{abc}+t_i^at_{jk}^{bc}-t_j^at_{ik}^{bc}+t_k^at_{ij}^{bc}-t_i^bt_{jk}^{ac}+t_j^bt_{ik}^{ac}-t_k^bt_{ij}^{ac}+t_i^ct_{jk}^{ab}-t_j^ct_{ik}^{ab}\\
&+t_k^ct_{ij}^{ab}+t_i^at_j^bt_k^c-t_i^at_j^ct_k^b-t_i^bt_j^at_k^c+t_i^bt_j^ct_k^a+t_i^ct_j^at_k^b-t_i^ct_j^bt_k^a\\
\ldots
\end{aligned}
\end{equation}
$$

由此可知，算符 $\hat{C}$ 中的轨道系数与CC振幅存在一一对应关系。因此，我们可以将式（17）写成由簇算符 $\hat{T}$ 表示的形式，即

$$
\begin{equation}
\ket{\Psi} = e^{\hat{T}}\ket{\Psi_0}
\end{equation}
$$

这就是精确波函数的耦合簇展开。从上述推导过程中可以看到，精确波函数的簇展开与 CI 展开实际上是等价的。将波函数展开成（17）的计算方法被称为**组态相互作用（Configure Interaction，CI）理论**，而将波函数展开成（23）的计算方法被称为**耦合簇（Coupled Cluster，CC）理论**。

### 3.2 CC能量方程与振幅方程

&emsp;&emsp;在上一小节，我们废了老大劲才从CI展开式（17）过渡到CC展开式（23），一些初学者可能会感到困惑：为什么要把原本的激发算符给转化成指数算符（23）呢，这种指数形式的算符处理起来与原来的激发算符相比更为简单吗？乍一看这一通操作好像把原来展开式变得更抽象更复杂了，但实际上，将算符转化成指数形式将会大大简化后续的推导。在这一小节中我们将会推导CC能量方程与振幅方程，在推导过程中我们就能体会到指数算符的便利之处。

&emsp;&emsp;在将精确波函数写成式（23）的形式后，因为参考态一般是给定的，我们只需做个 HF 单点计算就能得到 HF 基态能量以及基态轨道系数和密度矩阵，所以实际上我们只需要知道振幅的值就能够确定精确波函数的表达式，进而计算能量以及一些其他的我们感兴趣的物理量，从而确定体系的所有性质。

&emsp;&emsp;下面我们就来推导CC能量和振幅所满足的方程。将（23）代入定态薛定谔方程，可得

$$
\begin{equation}
\hat{H}e^{\hat{T}}\ket{\Psi_0} = Ee^{\hat{T}}\Psi_0
\end{equation}
$$

两边同时左乘 $e^{-\hat{T}}$，可得

$$
\begin{equation}
e^{-\hat{T}}\hat{H}e^{\hat{T}}\ket{\Psi_0}=\hat{H}_{cc}\ket{\Psi_0}=E\ket{\Psi_0}
\end{equation}
$$

其中  

$$
\begin{equation}
\hat{H}_{cc}=e^{-\hat{T}}\hat{H}e^{\hat{T}}
\end{equation}
$$

是哈密顿算符 $\hat{H}$ 的相似变换。在得到方程（25）后，两边左乘 $\bra{\Psi_0}$ ，根据组态函数的正交归一性，可以得到 CC 能量的表达式，即

$$
\begin{equation}
E = \bra{\Psi_0}\hat{H}_{cc}\ket{\Psi_0}
\end{equation}
$$

同理，两边同时左乘 $\bra{\Psi_i^a}$ 即可得到 $t_1$ 振幅所满足的方程，即

$$
\begin{equation}
0 = \bra{\Psi_i^a}\hat{H}_{cc}\ket{\Psi_0}
\end{equation}
$$

同理，两边同时左乘 $\bra{\Psi_{ij}^{ab}}$ 即可得到 $t_2$ 振幅所满足的方程，即

$$
\begin{equation}
0 = \bra{\Psi_{ij}^{ab}}\hat{H}_{cc}\ket{\Psi_0}
\end{equation}
$$

以此类推，最终我们可以得到任意阶振幅所满足的方程，即

$$
\begin{equation}
\begin{cases}
\bra{\Psi_i^a}\hat{H}_{cc}\ket{\Psi_0} & =0\\
\bra{\Psi_{ij}^{ab}}\hat{H}_{cc}\ket{\Psi_0} &=0\\
\bra{\Psi_{ijk}^{abc}}\hat{H}_{cc}\ket{\Psi_0} &=0\\
\qquad\cdots
\end{cases}
\end{equation}
$$

现在能量方程（27）和振幅方程（30）中只剩下 $\hat{H}_{cc}$ 未能显式表示。为了计算 $\hat{H}_{cc}$ ，我们需要引入一个公式——Baker–Campbell–Hausdorff 公式，其表达式如下（证明过程见附录A）

$$
\begin{equation}
e^{-A}Be^{A} = B+[B,A]+\frac{1}{2!}[[B,A],A]+\frac{1}{3!}[[[B,A],A],A]+\cdots
\end{equation}
$$

将上式代入方程（26）中，可得

$$
\begin{equation}
\hat{H}_{cc}=e^{-\hat{T}}\hat{H}e^{\hat{T}}= \hat{H}+[\hat{H},\hat{T}]+\frac{1}{2!}[[\hat{H},\hat{T}],\hat{T}]+\frac{1}{3!}[[[\hat{H},\hat{T}],\hat{T}],\hat{T}]+\cdots
\end{equation}
$$

这就是 $\hat{H}_{cc}$ 的表达式。将上式代入能量方程（27）和振幅方程（30）中，我们就能得到能量的具体表达式以及各阶振幅所满足的方程组。

### 3.3 讨论

&emsp;&emsp;通过上面的推导过程，我们就能很清楚地看到指数算符的优越性。如果我们直接将 CI 的激发算符（18）代入定态薛定谔方程，我们就会发现我们很难直接求得算符 $\hat{C}$ 的逆，所以在 CI 中，我们一般不讲精确波函数写成用激发算符作用在参考态上的形式，而是直接将式（2）代入定态薛定谔方程中，直接求解本征值方程，最后可以化简成一个对角化问题[^4]。而 CC 由于使用的是指数形式的激发算符，能够直接使用 BCH 公式 （31），不需要求出激发算符的逆就能得到 $\hat{H}_{cc}$ 的表达式，从而将问题给转化成求解振幅方程组（30）[^5]。

&emsp;&emsp;如果我们将 CI 展开式（2）和簇算符（19）展开至无穷阶，即使用Full CI或Full CC，那么 CI 和 CC 等价，这里上文已经证明过，但是在实际计算中，在大多数情况下，随着展开式阶数的增加，所需的计算时间会呈几何式增加，所以我们一般会对展开式进行截断，从而节省计算时间。对于 CI ，我们一般截断展开式（2），如 CISD 代表在展开式（2）中只保留基态、单激发组态和双激发组态，双激发以上的组态会被舍弃，即令

$$
\begin{equation}
\ket{\Psi} = \ket{\Psi_0}+\sum_{ia}c_i^a\ket{\Psi_i^a}+\sum_{iajb}c_{ij}^{ab}\ket{\Psi_{ij}^{ab}}
\end{equation}
$$

而在 CC 中，我们一般会截断簇算符的表达式（19），如 CCSD 代表在簇算符中只保留单激发和双激发部分，其他激发算符会被舍弃，即令

$$
\begin{equation}
\begin{aligned}
\hat{T} &= \hat{T}_1+\hat{T}_2\\
\ket{\Psi} &=e^{\hat{T}}\ket{\Psi_0}= e^{\hat{T}_1+\hat{T}_2}\ket{\Psi_0}
\end{aligned}
\end{equation}
$$

在这里，一些读者可能就要说了：在对簇算符进行截断后，式（27）（30）不还是包含 $\hat{H}_{cc}$ 这个无穷阶求和吗，这不还是需要进行大量计算吗？（32）看似是无穷阶求和，但实际上在计算（27）（30）等积分的时候会被自然截断，这是因为分子体系的哈密顿算符最多是二体算符，所以在实际计算的时候（27）（30）都会被会被自然截断至有限阶项。这一点在下一篇博文中将会进行详细讨论。

&emsp;&emsp;在截断的 CI 中，我们仅保留低阶的激发组态，把高阶激发组态全部舍弃，即假设高阶激发组态的组态系数远小于低阶激发组态的组态系数，这种截断方案其实并不妥当，因为这个假设并不一定成立，高阶组态系数未必一定会比低阶组态系数小，所以使用截断 CI计算得到的相关能精度十分有限。而截断的 CC 并不存在这个问题。由于在截断簇算符后，$\hat{H}_{cc}$ 中仍然会包含高阶激发态，所以 CC 的截断相较于 CI 的显得更为合理。由于相同阶数的截断CI和CC耗时相近看，因此在相同的计算时间下，CC得到的相关能一般会比 CI 的更加准确。

## 附录A. Baker–Campbell–Hausdorff 公式的简单证明

&emsp;&emsp;这里的证明步骤参考了[^6]。定义算符函数 $U(x)$

$$
\begin{equation}
U(x) = e^{-Ax}Be^{Ax}
\end{equation}
$$

其中 $x$ 为实数。上式两边对 $x$ 求导，可以求得 $U(x)$ 所满足的微分方程，即

$$
\begin{equation}
\frac{d}{dx}U=-Ae^{-Ax}Be^{Ax}+e^{-Ax}BAe^{Ax}
\end{equation}
$$

由于 算符 $A$ 与它的指数算符对易[^7]，上式可简化为

$$
\begin{equation}
\frac{d}{dx}U=-AU+UA=[U,A]
\end{equation}
$$

这就是算符函数 $U(x)$ 所满足的微分方程。两边同时积分，使用初值条件 $U(0)=B$，即可得到方程（37）的形式解，即

$$
\begin{equation}
U(x) = B+\int_0^x[U(t),A]dt
\end{equation}
$$

将 $U(x)$ 对 $x$ 进行泰勒展开，可得

$$
\begin{equation}
\begin{aligned}
\text{左边}&=U_0+U_1x+\frac{1}{2!}U_2x^2+\cdots\\
\text{右边}&=B+\int_0^xU(t)Adt-\int_0^xAU(t)dt\\
&=B+[U_0,A]x+\frac{1}{2!}[U_1,A]x^2+\cdots
\end{aligned}
\end{equation}
$$

对比两边系数，可得

$$
\begin{equation}
\begin{aligned}
U_0 &= B\\
U_1 &=[U_0,A]=[B,A]\\
U_2 &=[U_1,A]=[[B,A],A]\\
&\cdots
\end{aligned}
\end{equation}
$$

将系数代入（39）中，即可得算符函数 $U(x)$ 的级数表达式，即

$$
\begin{equation}
\begin{aligned}
U(x)&=U_0+U_1x+\frac{1}{2!}U_2x^2+\cdots=e^{-Ax}Be^{Ax} \\
&= B+[B,A]x+\frac{1}{2!}[[B,A],A]x^2+\frac{1}{3!}[[[B,A],A],A]x^3+\cdots
\end{aligned}
\end{equation}
$$

令 $x=1$，即可得到 BCH 公式（31），即

$$
\begin{equation}
\begin{aligned}
U(1)&=U_0+U_1+\frac{1}{2!}U_2+\cdots=e^{-A}Be^{A} \\
&= B+[B,A]+\frac{1}{2!}[[B,A],A]+\frac{1}{3!}[[[B,A],A],A]+\cdots
\end{aligned}
\end{equation}
$$

[^1]: <https://www.wikiwand.com/en/articles/Coupled_cluster>
[^2]: 这只是算符收缩的一种写法，更常见的写法是在算符的上方用形如中括号的折线连接待收缩的两个算符，但是markdown并不支持这个符号，所以这里使用另一种标记方法，即将两个上标上有相同数目的圆点的算符进行收缩。
[^3]: <https://zhuanlan.zhihu.com/p/675266512。>
除此之外徐光宪的量化下册p39也有wick定理的证明过程。
[^4]: szabo 书中有详细讨论。
[^5]: 该方程组并不是线性方程组，因为其中存在高阶项。实际上 CC 也可以将问题转化成类似与 CI 的本征值问题，但是由于构建的哈密顿矩阵并不对称，所以处理起来比较麻烦。一般程序做 CC 单点计算都是直接求解方程（30），但是 CC 的哈密顿矩阵还有其他用途，比如在做 CC 梯度计算的时候需要用到 CC 的哈密顿矩阵。这方面作者目前了解不多，这里就不过多讨论了。
[^6]: <https://zhuanlan.zhihu.com/p/618624483>
[^7]: 将指数算符进行泰勒展开后容易看出它们对易.
