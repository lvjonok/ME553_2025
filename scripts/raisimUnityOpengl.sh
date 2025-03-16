uname_os="$(uname -s)"
uname_machine="$(uname -m)"
case "${uname_os}" in
    Linux*)    ./third_party/raisimLib/raisimUnityOpengl/linux/raisimUnity.x86_64;;
    Darwin*)    echo "OpenGl Unity is not defined for Mac, please use raisimUnity.sh"; exit 1;;
    *)          echo "machine ${uname_os} is undefined"; exit 1;;
esac
