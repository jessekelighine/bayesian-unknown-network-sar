" format.vim

silent :RemoveTrailingSpaces

" remove first list
silent g/post\.mean/d

" format numbers
silent %s/-\?[0-9]\+\.[0-9]\+\>/$\0$/g

" format parameters
silent %s/beta\(\d\+\)/$\\bbeta_{\1}$/g
silent %s/lambda/$\\lambda$/g
silent %s/sigma2/$\\sigma^{2}$/g
silent %s/theta\(\d\+\)/$\\ttheta_{\1}$/g

" insert table separators
silent %s/\s\+/ \& /g
silent %s/$/ \\\\/

" align
silent execute "norm ggvipga*&\<CR>"
