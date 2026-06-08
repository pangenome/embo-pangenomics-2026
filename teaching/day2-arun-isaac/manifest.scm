(use-modules (guix download)
             (guix git-download)
             (guix packages)
             ((gnu packages bioinformatics)
              #:select (ccwl ravanan) #:prefix guix:)
             ((gnu packages graphviz) #:select (graphviz))
             ((gnu packages guile) #:select (guile-3.0-latest))
             ((gnu packages guile-xyz) #:select (guile-run64)))

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

(define-public ccwl
  (package
    (inherit guix:ccwl)
    (name "ccwl")
    (version "0.5.0")
    (source
     (origin
       (method url-fetch)
       (uri (string-append "https://ccwl.systemreboot.net/releases/ccwl-"
                           version ".tar.lz"))
       (sha256
        (base32
         "1pna61cqhhg69qpa682x911x85nlin1wpp29naqzq3r0s951qxi6"))))
    (native-inputs
     (modify-inputs (package-native-inputs guix:ccwl)
       (prepend guile-run64)))))

(packages->manifest
 (list ccwl graphviz ravanan))
