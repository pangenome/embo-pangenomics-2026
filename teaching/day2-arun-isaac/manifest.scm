(use-modules (guix git-download)
             (guix packages)
             ((gnu packages bioinformatics) #:select (ccwl))
             ((gnu packages bioinformatics)
              #:select (ravanan) #:prefix guix:)
             ((gnu packages graphviz) #:select (graphviz))
             ((gnu packages guile) #:select (guile-3.0-latest)))

(define-public ravanan
  (let ((commit "1f5000ad6ff98278bf638ff0176ddfd5bf8933cf"))
    (package
      (inherit guix:ravanan)
      (name "ravanan")
      (version "0.2.0")
      (source (origin
                (method git-fetch)
                (uri (git-reference
                      (url "https://git.systemreboot.net/ravanan")
                      (commit commit)))
                (file-name (git-file-name name version))
                (sha256
                 (base32
                  "1dq3jn8z78x287krfvi6hbj066pmz8mk2qamjfsind0bxxjk3qzd"))))
      (inputs
       (modify-inputs (package-inputs guix:ravanan)
                      (replace "guile" guile-3.0-latest))))))

(packages->manifest
 (list ccwl graphviz ravanan))
