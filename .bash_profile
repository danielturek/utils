
# to permanently mute OSX startup chime, run this once from shell:
# $ cd ~/github/utils/nobootsound
# $ sudo sh install.sh

# download and install f.lux program, for night shift softening of colours.

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

alias e=emacs
alias ex=exit
alias o='open -a "Google Chrome"'
alias rmm=~/github/utils/move_to_trash.sh
alias ls='ls -F'
alias backup=~/github/utils/backup.sh
alias cleanup=~/github/utils/cleanup.sh
alias build=~/github/utils/build.sh
alias gand='ssh gandalf.berkeley.edu'
alias hpcc='ssh dbt1@hpcc.williams.edu'
alias status='git status'
alias st='git status'
alias dif='git diff'
alias pull='git pull'
alias push='git push'
alias checkout='git checkout'
alias branch='git branch'
alias commit=~/github/utils/commit.sh


set -o noclobber

