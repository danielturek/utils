# to permanently mute OSX startup chime, run this once from shell:
# $ cd ~/github/utils/nobootsound
# $ sudo sh install.sh

# download and install f.lux program, for night shift softening of colours.

# download and install Jumpcut program, for better copy/paste options

# install homebrew:
# $ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# install gsed (GNU sed):
# $ brew install gnu-sed

# stop creation of .DS_Store files by Finder everywhere
# $ defaults write com.apple.desktopservices DSDontWriteNetworkStores true

# to turn on "3-finger drag", under
# System Prefereces, Accessibility, Mouse & Trackpad, Trackpad Options, Enable dragging

# need to download & install to get jags to run on El Capitan:
# JAGS-Mavericks-3.4.0.dmg
# then do usual install.packages('rjags')

# download and install MacTeX (which is huge), then have to update paths in:
# TeXShop, LaTeXiT, BibDesk, TeX Live Utility, to /Library/TeX/texbin
# see: https://tug.org/mactex/UpdatingForElCapitan.pdf

# download and install Anaconda Python version 3.7, from:
# (the "command line" version of the installer works better - use that one)
# https://www.anaconda.com/distribution/
# after the Anaconda insallation, run from terminal (to remove "(base)" from bash prompt):
# $ conda config --set changeps1 False

# to get polymode (enter into R style-editing in Rmarkdown files),
# accessed by mmy own shortcut: C-c C-m, in Aquamacs, run:
# M-x package-install [ENTER] polymode      [ENTER]
# M-x package-install [ENTER] poly-markdown [ENTER]
# M-x package-install [ENTER] poly-R        [ENTER]

# Macbook FIXES
# 
# Apple Technical Support (Education): 1-800-800-2775, option 3
# 5/1/2018 case # 100523958608
#
# force immediate shutdown: shit+control+option+power
#
# built-in keyboard and trackpad not working whatsoever:
# (1) with computer on and screen open,
# "pinch" the middle-right-side of the computer,
# to reconnect loose pins....
# (2) shut down computer,
# press (left side) shift+control+option+power for a while,
# then release them at the same time,
# then turn computer on using power key.
# (3) shut down computer,
# hold down (left side) command+option+P+R, then
# turn on computer using power key,
# continue holding them down until computer starts up,
# then release.

# to turn off "red notification circle" appearing on
# the Systems Preferences icon in the Dock
# (which is reminding you to upgrade the OS):
# $ defaults write com.apple.systempreferences AttentionPrefBundleIDs 0
# $ killall Dock

# in PS1 variable to define prompt:
# see: http://www.funtoo.org/Prompt_Magic
# "\a" makes a little bell noise
# "\w" current working directory
# "\u" username
# "\h" hostname (\H for long hostname)
# "\d" date
# "\@" time
# "\e[...m" is the enclosing sequence for colors, multiple values in ... separated by ";"
PS1="\[\e[36m\]\u@\h on \d at \@\n\w\[\e[32m\] \$(parse_git_branch)\n\[\e[0m\]$ "

parse_git_branch() {
    git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/^* \(.*\)/(branch: \1)/'
}

function g() {
  cd ~/github/$1
}

alias e='emacs -nw'
alias ex=exit
alias o='open -a "Google Chrome"'
alias textbook=~/github/utils/textbook.sh
alias solutions='open ~/github/courses/klingenberg_textbook/v5_textbook_solutions/v5_solutions.pdf'
alias lecture=~/github/utils/lecture.sh
alias rmm=~/github/utils/move_to_trash.sh
alias ls='ls -F'
alias cv='open ~/github/private/CV/Turek\ -\ CV.tex'
alias cvcommit='dir=$(pwd); cd ~/github/private; commit; cd $dir'
alias backup=~/github/utils/backup.sh
alias cleanup=~/github/utils/cleanup.sh
alias build=~/github/utils/build.sh
alias buildMatt=~/github/utils/build_matt.sh
alias Rrun=~/github/utils/run_R_script.sh
alias gand='ssh dturek@gandalf.berkeley.edu'
alias hpcc='ssh dbt1@hpcc.williams.edu'
alias status='git status'
alias st='git status'
alias dif=~/github/utils/my_github_dif.sh
alias pull='git pull'
alias push='git push'
alias checkout='git checkout'
alias branch='git branch'
alias log='git log'
alias commit=~/github/utils/commit.sh
alias makesite=~/github/utils/make_site.sh
## alias s101='cd ~/github/courses/stat101'
## alias s161='cd ~/github/courses/stat161'
## alias s201='cd ~/github/courses/stat201'
## alias s202='cd ~/github/courses/stat202'
## alias s319='cd ~/github/courses/stat319'
## alias s359='cd ~/github/courses/stat359'
## alias s360='cd ~/github/courses/stat360'
## alias s365='cd ~/github/courses/stat365'
## alias s465='cd ~/github/courses/stat465'
##
alias nim='cd ~/github/nimble/nimble/packages/nimble/R'
##alias m='cd ~/github/courses/math336'
##alias m1='cd ~/github/courses/math182'
##alias m3='cd ~/github/courses/math335'
alias m3='cd ~/github/courses/math371'
alias ds='cd ~/github/ds201/ds201'
alias site='cd ~/github/ds201/site'
alias slides='cd ~/github/ds201/ds201/slides'
alias starters='cd ~/github/ds201/starters'
##alias edit=~/github/utils/edit.sh   doesn't work
##alias stop=~/github/utils/stop.sh   doesn't work
alias makeassignment=~/github/utils/make_assignment.sh

set -o noclobber

PATH="$PATH:/opt/homebrew/bin/"

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
#__conda_setup="$('/Users/dturek/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
#if [ $? -eq 0 ]; then
#    eval "$__conda_setup"
#else
#    if [ -f "/Users/dturek/opt/miniconda3/etc/profile.d/conda.sh" ]; then
#        . "/Users/dturek/opt/miniconda3/etc/profile.d/conda.sh"
#    else
#        export PATH="/Users/dturek/opt/miniconda3/bin:$PATH"
#    fi
#fi
#unset __conda_setup
# <<< conda initialize <<<

