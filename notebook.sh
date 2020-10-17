##### Customized
# run this script from any directory 
# it will load anaconda from the specified directory
ROOT_PATH="$(echo) $HOME/softwares"
echo $ROOT_PATH
command="$ROOT_PATH/anaconda3/bin/conda"

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$($command 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "$ROOT_PATH/anaconda3/etc/profile.d/conda.sh" ]; then
        . "$ROOT_PATH/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="$ROOT_PATH/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


path_to_notebook=$(pwd)
echo $path_to_notebook
jupyter notebook $path_to_notebook

