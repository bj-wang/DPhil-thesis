%% Section 0: Declaration
%% 第〇节：声明
% 
% LaTeX class template of the DPhil thesis, Bingjun Wang, Department of Physics, University of Oxford, May 2021
% 博士论文 LaTeX 文类模板，王炳俊，牛津大学物理系，2021 年 5 月
% 
% Main references: (1) John McManigle. ociamthesis - Oxford thesis template. https://www.oxfordechoes.com/oxford-thesis-template/.
% Main references: (2) Wei Hu. A Comprehensive Handbook of LaTeX2e (2nd edn., in Chinese). Tsinghua University Press, 2013.
% Main references: (3) Wei Hu. A Handbook of Classes and Packages in LaTeX2e (in Chinese). Tsinghua University Press, 2017.
% 主要参考文献：（1）John McManigle. ociamthesis - Oxford thesis template. https://www.oxfordechoes.com/oxford-thesis-template/.
% 主要参考文献：（2）胡伟. LaTeX2e 完全学习手册（第二版）. 清华大学出版社, 2013.
% 主要参考文献：（3）胡伟. LaTeX2e 文类和宏包学习手册. 清华大学出版社, 2017.
% 
% This file is licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License: https://creativecommons.org/licenses/by-nc-nd/4.0/
% 本文件受到署名-非商业性使用-禁止演绎 4.0 国际知识共享协议的保护：https://creativecommons.org/licenses/by-nc-nd/4.0/deed.zh
% 
% Please contact the author (bingjun.wang@physics.ox.ac.uk or bingjun.wang1995@outlook.com) if you have any questions or comments.
% 如有问题或意见，请联系作者：bingjun.wang@physics.ox.ac.uk 或 bingjun.wang1995@outlook.com
% 
% 
%% Section 1: Indentification
%% 第一节：版本设定
% 
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{BW-DPhil-thesis-class}[2021/05/14 BW DPhil thesis class]
\LoadClass[12pt,a4paper]{book}
% 
% 
%% Section 2: Option Declaration
%% 第二节：文类选项声明
% 
% Define a Boolean register to decide whether this is the printed version. Default value is false.
% 定义一个布尔寄存器来判断是否为打印版本，默认值为否
\newif\ifprintedver
\printedverfalse
% Apply the printedver option to the class.
% 将 printedver 选项应用到文类中
\DeclareOption{printedver}{\printedvertrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ProcessOptions\relax
% 
% 
%% Section 3: Package Loading
%% 第三节：导入宏包
% 
% Section 3.1: Page layout    % 第 3.1 节：版面设计宏包
% \RequirePackage{geometry}    % set margins, has to be loaded after the hyperref package    % 设置页边距，需在 hyperref 宏包之后导入
\RequirePackage{setspace}    % set line spacing    % 设置行距
\RequirePackage{fancyhdr}    % set headers and footers    % 设置页眉页脚
% 
% Section 3.2: Environments    % 第 3.2 节：环境宏包
\RequirePackage{graphicx}    % figure insertion    % 插图
\RequirePackage{booktabs}    % three-line tables    % 三线表
\RequirePackage{threeparttable}    % footnotes for tables    % 表格脚注
\RequirePackage{makecell}    % forced new lines in tables    % 表格中的强制换行
\RequirePackage{paralist}    % a better list environment    % 更好的列表环境
\RequirePackage[titles]{tocloft}    % fine tune the list of figures and tables    % 对图表目录进行微调
\RequirePackage{datatool}    % sort the list of abbreviations and symbols    % 对缩写和符号表进行排序
\RequirePackage[style=nature]{biblatex}    % generate the bibliography    % 产生参考文献
% 
% Section 3.3: Math and chemistry    % 第 3.3 节：数学和化学宏包
\RequirePackage{amsmath}    % bacic math    % 基本数学宏包
\RequirePackage{amssymb}    % bacic math    % 基本数学宏包
\RequirePackage{amstext}    % text in the math mode    % 数学模式中的文本
\RequirePackage{siunitx}    % numbers and units in the main text    % 正文中的数字和单位
\RequirePackage{esvect}    % better arrows for vectors    % 更好的向量箭头
\RequirePackage{unicode-math}    % upright circular constant and partial derivative symbols    % 直立的圆周率和偏微分符号
\setmathfont{latinmodern-math.otf}    % set the math font; has to be placed right after the loading of the unicode-math package    % 设置数学字体，必须紧跟在导入 unicode-math 宏包之后
\RequirePackage{chemformula}    % chemical formulas    % 化学式
\RequirePackage{chemfig}    % chemical structures    % 化学结构式
\RequirePackage{chemmacros}    % long chemical names    % 长化学名
% 
% Section 3.4: Titles and captions    % 第 3.4 节：标题宏包
\RequirePackage{titlesec}    % chapter and section titles    % 章节标题
\RequirePackage{caption}    % figure and table captions    % 图表标题
% 
% Section 3.5: Miscellaneous    % 第 3.5 节：杂项宏包
\RequirePackage{xcolor}    % colours    % 颜色宏包
\RequirePackage[EU1]{fontenc}    % accented characters    % 重音字符
\RequirePackage{fontspec}    % choose font    % 选择字体
\RequirePackage{microtype}    % fine-tuning the letter intervals    % 字距微调
\RequirePackage{calc}    % calculation in the \setlength command    % \setlength 命令中的计算
\RequirePackage{xfakebold}    % generate fake bold words    % 伪加粗字体
\RequirePackage{xparse}    % provides a \NewDocumentEnvironment command that allows auguments to be passed to the end part of a new envrionment    % 提供了 \NewDocumentEnvironment 命令，使得新环境的结束定义可以使用参数 
% 
% Section 3.6: Citation    % 第 3.6 节：引用宏包
\RequirePackage{hyperref}    % cross-referencing and hyperlinks    % 交叉引用和超链接
\RequirePackage[capitalise,nosort]{cleveref}    % better cross-referencing    % 更好的交叉引用
% 
% Section 3.7: Conflicts - have to be placed here   % 第 3.7 节：冲突宏包，只能置于此
\RequirePackage{geometry}
% 
% 
%% Section 4: Detailed Page Settings
%% 第四节：详细页面设置
% 
% Section 4.1: Page layout (depending on being printed or not)    % 第 4.1 节：版面，与是否打印有关
\ifprintedver
	\geometry{includehead,hmargin={3cm,2cm},vmargin={2.5cm,2.5cm}}
\else
	\geometry{includehead,hmargin={2.5cm,2.5cm},vmargin={2.5cm,2.5cm}}
\fi
% 
% Section 4.2: Line spacing    % 第 4.2 节：行距
\doublespacing
% 
% Section 4.3: Headers and footers    % 第 4.3 节：页眉页脚
\renewcommand\chaptermark[1]{\markboth{\thechapter. \ #1}{}}    % set the format of chapter titles in headers    % 设置页眉区章标题的格式
\pagestyle{fancy}
\fancyhf{}
\renewcommand{\headrulewidth}{0pt}    % no head rule    % 无页眉线
\fancyheadoffset[EL,OR]{3em}    % set header offset to accommodate page numbers    % 页眉区向外凸出以排版页码
\fancyhead[EL]{\makebox[2em][r]{\textsl{\thepage}}\quad|\quad\nouppercase{\textsl{\leftmark}}}    % even page left headers    % 偶数页左页眉
\fancyhead[OR]{\nouppercase{\textsl{\rightmark}}\quad|\quad\makebox[2em][l]{\textsl{\thepage}}}    % odd page right headers    % 奇数页右页眉
% Re-define the plain style    % 重定义 plain 版式
\fancypagestyle{plain}{%
	\fancyhf{}%
	\fancyhead[EL]{\makebox[2em][r]{\textsl{\thepage}}\quad}%
	\fancyhead[OR]{\quad\makebox[2em][l]{\textsl{\thepage}}}%
}
% Re-define the \cleardoublepage command for the (possible) blank page after a chapter    % 重定义 \cleardoublepage 清双页命令，应用于可能存在的章末空白页
\renewcommand*{\cleardoublepage}{%
	\clearpage%
		\if@twoside%
			\ifodd\c@page%
		\else%
			\hbox{}%
			\thispagestyle{plain}%
			\newpage%
			\if@twocolumn%
				\hbox{}%
				\newpage%
			\fi%
		\fi%
	\fi%
}
\addtolength{\topmargin}{-2.5pt}    % reduce the top margin    % 减少顶部空白
\setlength{\headheight}{14.5pt}    % set the header height    % 设置页眉区高度
% Enable short section title in the header    % 页眉中的短节标题
\newcommand{\presecmarkBW}[1]{%
	\addtocounter{section}{1}%
	\sectionmark{#1}%
	\addtocounter{section}{-1}%
}
% 
% Section 4.4: Chapter titles    % 第 4.4 节：章标题
% For numbered chapters    % 有编号章
\titleformat{name=\chapter}[display]%
	{\bf\Large}%
	{\filleft\sc\LARGE{\chaptertitlename} \fontsize{50}{60}\usefont{OT1}{uaq}{m}{n}\selectfont\thechapter}%
	{1ex}%
	{\titlerule[1pt]\vspace{0.5ex}\centering}%
	[\vspace{1ex}{\titlerule[1pt]}]%
% For non-numbered chapters    % 无编号章
\titleformat{name=\chapter,numberless}[display]
	{\sc\LARGE}%
	{}%
	{0pt}%
	{\filleft}%
	[\vspace{1ex}]%
% Surrounding spacing of the chapter title    % 章标题四周间距
\titlespacing{\chapter}{0pt}{*-4}{*3}
% 
% Section 4.5: Section titles    % 第 4.5 节：节标题
\setcounter{secnumdepth}{3}    % numbering the subsubsection title    % 给出小小节标题的编号
\titleformat*{\section}{\large\bf}
\titleformat*{\subsection}{\normalsize\bf}
\titleformat*{\subsubsection}{\normalsize\sl}
\titleformat{\paragraph}[hang]{\normalsize\sl}{}{0pt}{}
\titlespacing{\paragraph}{2em}{*2}{*2}
% 
% Section 4.6: Captions    % 第 4.6 节：图表标题
\DeclareCaptionLabelSeparator{vertline}{ | }    % define the new label separator    % 定义新的标签分隔符
% Set caption styles, including label separator, font, line spacing, and surrounding spacing    % 设置图表标题格式，包括标签分隔符、字体属性、行距、上下文间距
\captionsetup{%
	format=plain,%
	labelsep=vertline,%
	labelfont=bf,%
	font={onehalfspacing,small},%
	singlelinecheck=false,%
	belowskip=5pt%
}
% 
% Section 4.7: Fine tuning    % 第 4.7 节：微调
\setmainfont{Latin Modern Roman}[SmallCapsFont={* Caps},SlantedFont={* Slanted}]    % select fonts    % 选择字体
\setlength{\tabcolsep}{1em}    % half of the separation between table columns    % 表格各列的半间距
\renewcommand{\arraystretch}{1.5}    % table stretching coefficient    % 表格拉伸系数
\setlength{\intextsep}{15pt plus 2pt minus 2pt}    % separation between floats and the main text    % 浮动体与正文的间距
\pretolerance=1000    % tolerance of badness (w/o hyphenation), default = 100    % 无断词时的劣质参数容许度，默认为 100
\tolerance=2000    % tolerance of badness (w/ hyphenation), default = 200    % 有断词时的劣质参数容许度，默认为 200
\hyphenpenalty=1000    % the ability to inhibit hyphenation, default = 50    % 抑制断词的能力，默认为 50
\lefthyphenmin=3    % when hyphenation, the minimum number of letters on the left of the hyphen, default = 2    % 断词时连字符左侧的最小字母数，默认为 2
\righthyphenmin=4    % when hyphenation, the minimum number of letters on the right of the hyphen, default = 3    % 断词时连字符右侧的最小字母数，默认为 3
\widowpenalty=500    % the ability to inhibit widow lines, default = 150    % 抑制页首孤行的能力，默认为 150
\clubpenalty=500    % the ability to inhibit orphan lines, default = 150    % 抑制页尾孤行的能力，默认为 150
%
% Section 4.8: Package options    % 第 4.8 节：宏包选项
% hyperref package: link colours, and numbers in the pdf bookmark    % hyperref 宏包：链接颜色与 pdf 文件中的书签序号
\hypersetup{%
	colorlinks,%
	linkcolor=blue,%
	anchorcolor=black,%
	citecolor=blue,%
	urlcolor=blue,%
	bookmarksnumbered%
}
% cref package: the conjunction between the last two items of a list (add the Oxford comma)    % cref 宏包：引用项列表的最后两项之间的连接词（增加了连接逗号）
\newcommand{\creflastconjunction}{, and }
% siunitx package: units and conjuctions in lists, products, and ranges; show explicit plus sign; do not show ``1\times'' in the scientific notation of numbers; show uncertainty as \pm    % siunitx 宏包：设置 list、product 和 range 中的单位和连接词，按输入方式显示正号，不显示科学计数法中的``1\times''，以 \pm 符号的方式显示不确定度
\sisetup{%
	range-units=single,%
	range-phrase={--},%
	list-units=single,%
	list-final-separator={, and },%
	product-units=single,%
	retain-explicit-plus=true,%
	print-unity-mantissa=false,%
	separate-uncertainty=true,%
	separate-uncertainty-units=single%
}
% chemfig package: atom separation, bond width, and better-looking bond-joint-points    % chemfig 宏包：原子间距，化学键粗细，化学键连接处的更好显示
\setchemfig{%
	atom sep=1.5em,%
	bond style={line width=0.8pt},%
	bond join=true%
}
% chemmacros package: cancel all spacing caused by hyphens and hyphenation indicators (|) in the nomenclature module    % chemmacros 宏包：取消 nomenclature 模块中所有由连字符和断词指示符（|）引起的空白
\chemsetup[nomenclature]{%
	hyphen-pre-space=0pt,%
	hyphen-post-space=0pt,%
	break-space=0pt%
}
% 
% 
%% Section 5: The Cover Page
%% 第五节：首页设置
% 
% Define the information that the author needs to input
% 定义作者需要提供的信息
\newcommand{\college}[1]{\gdef\@college{#1}}
\newcommand{\termyear}[1]{\gdef\@termyear{#1}}
\newcommand{\Oxfordlogopath}[1]{\gdef\@Oxfordlogopath{#1}}
% Define a \maketitleBW command to generate the cover page
% 定义一个 \maketitleBW 命令以生成首页
\newcommand{\maketitleBW}{\pdfbookmark[0]{Cover}{coverpage}%
	\begin{titlepage}%
		\begin{center}%
			\begin{onehalfspace}%
				{\setBold[0.6]\Huge\@title\unsetBold \par}%
				\vspace*{80pt}%
				{\includegraphics[width=0.25\linewidth]{\@Oxfordlogopath} \par}%
				\vspace*{80pt}%
			\end{onehalfspace}%
			\begin{spacing}{1}%
				{\Large\@author \par}%
				\vspace*{8pt}%
				{\large\@college \par}%
				\vspace*{8pt}%
				{\large University of Oxford \par}%
				\vspace*{60pt}%
				{\large A thesis submitted for the degree of \par}%
				\vspace*{8pt}%
				{\large\it Doctor of Philosophy \par}%
				\vspace*{30pt}%
				{\large\@termyear \par}%
			\end{spacing}%
		\end{center}\par
	\end{titlepage}%
	\global\let\maketitleBW\relax%
	\global\let\@title\@empty%
	\global\let\@Oxfordlogopath\@empty%
	\global\let\@author\@empty%
	\global\let\@college\@empty%
	\global\let\@termyear\@empty%
	\global\let\title\relax%
	\global\let\Oxfordlogopath\relax%
	\global\let\author\relax%
	\global\let\college\relax%
	\global\let\termyear\relax%
}
% 
% 
%% Section 6: The Abstract and Acknowledgement Environment
%% 第六节：摘要和致谢环境
% 
% Contains one argument to indicate the chapter title
% 含一个参数，章标题
\newenvironment{absackBW}[1]{%
	\thispagestyle{plain}%
	\phantomsection%
	\chapter*{#1}%
	\addcontentsline{toc}{chapter}{#1}%
}%
{}
% 
% 
%% Section 7: Contents
%% 第七节：目录设置
% 
\setcounter{tocdepth}{3}    % show subsection titles in contents    % 在目录中显示小节标题
\renewcommand{\@dotsep}{2}    % dot separation, in mu (=1/18 em)    % 指引圆点的间距，单位为 mu（=1/18 em）
% Combine the list of figures and the list of tables (to list of figures and tables)    % 将图片目录和表格目录合并为图表目录
\renewcommand{\ext@table}{lof}
\renewcommand{\listfigurename}{List of Figures and Tables}
% Show ``Figure'' and ``Table'' in the list of figures and tables,set the item indent    % 在图表目录中显示``Figure''和``Table''字样，设置条目缩进
\newlength{\loftlenBW}
\renewcommand{\cftfigpresnum}{\figurename\ }
\renewcommand{\cftfigaftersnum}{.}
\settowidth{\loftlenBW}{\cftfigpresnum\cftfigaftersnum}
\addtolength{\cftfignumwidth}{\loftlenBW}
\renewcommand{\cfttabpresnum}{\tablename\ }
\renewcommand{\cfttabaftersnum}{.}
\addtolength{\cfttabnumwidth}{\loftlenBW}
% 
% 
%% Section 8: The List of Abbreviations and Symbols Environment
%% 第八节：缩写和符号表环境
% 
% Define a \sortitemBW command to receive the input of each abbreviation or symbol item. It contains three arguments, i.e. the sorting label, the short version, and the long version.
% 定义一个 \sortitemBW 命令以接受每一条缩写或符号的输入。含三个参数，分别为排序标签、短形式和长形式
\newcommand{\sortitemBW}[3]{%
	\DTLnewrow{listBW}%
	\DTLnewdbentry{listBW}{label}{#1}%
	\DTLnewdbentry{listBW}{short}{#2}%
	\DTLnewdbentry{listBW}{long}{#3}%
}
% Define an abbsymBW environment to generate the list of abbreviations and symbols. It contains two arguments, i.e. the longest ``short'' version of the item, and the separation between the label and the item.
% 定义一个 abbsymBW 环境来生成缩写和符号表，包含两个参数，第一个为最长的短形式条目，第二个为标号与条目的间距
\NewDocumentEnvironment{abbsymBW}{mm}{%
	\DTLifdbexists{listBW}{\DTLcleardb{listBW}}{\DTLnewdb{listBW}}%
}%
{%
	\DTLsort*{label}{listBW}%
	\begin{list}{}{%
			\renewcommand\makelabel[1]{##1\hfil}%
			\settowidth\labelwidth{\makelabel{#1}}%
			\itemsep=-4.5pt\labelsep=#2%
			\setlength\leftmargin{\labelwidth+\labelsep}%
		}%
		\DTLforeach*{listBW}{\theShort=short,\theLong=long}{\item[\theShort] \theLong}%
	\end{list}
}
% Note: the sorting order is number - mixed Latin letters and Greek letters (case insensitive). For example, 1 - 2 - 3 - A and a- \alpha - B and b - \beta - ... For input details, see sorting-notes/sorting-notes.pdf
% 注：排序的顺序为：数字 - 不分大小写的拉丁和希腊字母。例如，1 - 2 - 3 - A和a - \alpha - B和b - \beta - ... 详情可见sorting-notes/sorting-notes.pdf
% 
% 
%% Section 9: Bibliography
%% 第九节：参考文献设置
% 
% Settings: show year only, show doi, enable back-referencing (do not compressing page numbers), show yyyy-mm-dd url date
% 设置：只显示年份，显示 doi，开启反向引用（不压缩页码列表），显示网页访问日期为年-月-日
\ExecuteBibliographyOptions{%
	date=year,%
	doi=true,%
	backref=true,%
	backrefstyle=none,%
	urldate=iso%
}
\DefineBibliographyStrings{english}{backrefpage={p.},backrefpages={pp.}}    % set the singular and the plural forms of the words before page numbers in back-referencing    % 设置反向引用中页码前的指示词，包括单数和复数两种情况
\setlength{\bibitemsep}{0pt plus 1pt}    % set vertical spacing between items    % 设置条目间的垂直间距
% 
% 
%% Section 10: Customised Macros for Convenience
%% 第十节：为输入方便而定义的宏
% 
% Section 10.1: SI units    % 第 10.1 节：SI 单位
\DeclareSIUnit{\muunit}{\square\cm\per\V\per\s}
\DeclareSIUnit{\mlBW}{ml}
\DeclareSIUnit{\angstrom}{\textup{\AA}}
\DeclareSIUnit{\rpmBW}{rpm}
% 
% Section 10.2: IUPAC nomenclature    % 第 10.2 节：IUPAC 命名
\NewChemIUPAC{\altBW}{\textit{alt}}
\NewChemIUPAC{\coBW}{\textit{co}}
\NewChemIUPAC{\secBW}{\textit{sec}}
\NewChemIUPAC{\bBW}{\textit{b}}
\NewChemIUPAC{\cBW}{\textit{c}}
\NewChemIUPAC{\fBW}{\textit{f}}
\NewChemIUPAC{\nBW}{\textit{n}}
% 
% Section 10.3: Others    % 第 10.3 节：其他
\newcommand{\andor}{and\slash or }
\newcommand{\symBW}[2]{\textit{#1}\textsubscript{#2}}    % physical quantities    % 物理量
% \newcommand{\TBDBW}{{\color{red} \textbf{TO BE DONE}}}    % notification for TBD during drafting, should be useless in the end    % 草拟时的待完成提示，最终应不需使用