\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{rhitbio}[2021 Bio Thesis Class]
\newcommand{\headlinecolor}{\normalcolor}
\RequirePackage{xcolor}
\definecolor{slcolor}{HTML}{882B21}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass{article}
\usepackage{tocloft}
\setlength{\cfttabindent}{0in}   
\setlength{\cftfigindent}{0in}   
\usepackage{ragged2e}
\setlength{\RaggedRightParindent}{\parindent}
\RaggedRight
\usepackage[font=doublespacing]{caption}
\usepackage{lipsum}  
\usepackage{enumitem}
\usepackage{tabto}
\usepackage{indentfirst}
\usepackage{datetime}
\usepackage{setspace}

\newcommand{\names}[3]{
\def\@authorfirst{#1}
    \def\@authormid{#2}
    \def\@authorlast{#3}
    \author{\@authorfirst\ \@authormid\ \@authorlast}
}

\newcommand{\submitdate}[2]{
    \def\@datemonth{#1}
    \def\@dateyear{#2}
}

\newcommand{\advisor}[1]{
    \def\@advisor{#1}
}

\newcommand{\degree}[2]{
    \def\@degree{#1}
    \def\@degreeshort{#2}
}

\newcommand{\subtitle}[1]{
    \def\@subtitle{#1}
}

\newcommand{\blankpage}{
    \newpage
    \thispagestyle{empty}
    \mbox{}
    \newpage
}

\renewcommand{\maketitle}{%
    \doublespacing
    \begin{center}
    \vspace*{\fill}
    \fontshape{n}\fontseries{m}
    \textbf{\@title:\\
    \@subtitle}\\
    A Thesis\\
    Submitted to the Faculty\\
    of\\
    Rose-Hulman Institute of Technology\\
    by\\
    \@author\\
    In Partial Fulfillment of the Requirements for the Degree of\\
    Bachelor of Science in \@degree\\
    \@datemonth\ \@dateyear
    \vspace*{\fill}\\
    \copyright \@dateyear\ \@author
    \end{center}
    \newpage
}

\newcommand{\abstractpreamble}{
    \noindent
    \@authorlast, \@authorfirst\ \@authormid\\
    B.S.\@degreeshort.
    Rose-Hulman Institute of Technology\\
    \@datemonth\ \@dateyear\\
    \@title: \@subtitle \\
    Thesis Advisor: \@advisor
    \newline 
}

\renewenvironment{abstract}
 {
  \begin{center}
  \bfseries \abstractname\vspace{-.5em}\vspace{0pt}
  \end{center}
  \relax
}

\newcommand{\setcommentspace}{
    \setlength{\textwidth }{5in}
    \setlength{\marginparwidth }{5cm}
}

\renewcommand{\normalsize}{\fontsize{9}{10}\selectfont}
\setlength{\textwidth}{17.5cm}
\setlength{\textheight}{22cm}