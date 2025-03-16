uname_os="$(uname -s)"
uname_machine="$(uname -m)"
case "${uname_os}" in
    Linux*)     ./third_party/raisimLib/raisimUnityOpengl/linux/raisimUnity.x86_64;;
    Darwin*)    
                if [ "$uname_machine" = "x86_64" ]; then
                    ./third_party/raisimLib/raisimUnity/mac/Contents/MacOS/raisimUnity;
                elif [ "$uname_machine" = "arm64" ]; then
                    ./third_party/raisimLib/raisimUnity/m1/RaiSimUnity.app/Contents/MacOS/RaiSimUnity;
                fi;;
    *)          echo "machine ${uname_os} is undefined"; exit 1;;
esac
