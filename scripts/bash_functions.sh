# Functions 
# Function to reverse complement a DNA sequence
reverse_complement() {
    local sequence="$1"
    local complement=""
    local len=${#sequence}

    for ((i = len - 1; i >= 0; i--)); do
        char="${sequence:$i:1}"
        case "$char" in
            [Aa]) complement+="T" ;;
            [Tt]) complement+="A" ;;
            [Cc]) complement+="G" ;;
            [Gg]) complement+="C" ;;
            [Rr]) complement+="Y" ;;
            [Yy]) complement+="R" ;;
            [Ss]) complement+="S" ;;
            [Ww]) complement+="W" ;;
            [Kk]) complement+="M" ;;
            [Mm]) complement+="K" ;;
            [Bb]) complement+="V" ;;
            [Dd]) complement+="H" ;;
            [Hh]) complement+="D" ;;
            [Vv]) complement+="B" ;;
            *) complement+="$char" ;;
        esac
    done

    echo "$complement"
}
# End of functions