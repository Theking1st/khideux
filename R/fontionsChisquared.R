#' Test du khi-deux pour les tableaux de contingence et les vecteurs de fréquences
#'
#' Cette fonction effectue un test du khi-deux pour les tableaux de contingence et les vecteurs de fréquences.
#'
#' @param x Un tableau de contingence ou un vecteur de fréquences observées.
#' @param y Un vecteur de fréquences observées (facultatif).
#' @param correct Logique, si TRUE, applique la correction de continuité de Yates pour les tableaux 2x2.
#' @param p Un vecteur de probabilités théoriques pour les fréquences observées.
#' @param rescale.p Logique, si TRUE, les probabilités sont mises à l'échelle pour totaliser 1.
#' @param simulate.p.value Logique, si TRUE, utilise une simulation pour calculer la valeur p.
#' @param B Nombre de réplicats pour la simulation de la valeur p.
#' @return Un objet de classe "htest" contenant les résultats du test du khi-deux.
#' @examples
#' # Exemple avec un tableau de contingence
#' dataE <- matrix(c(50, 100, 120, 80, 130, 20), nrow = 3, byrow = TRUE,
#'                dimnames = list(c("Primaire", "Secondaire", "Supérieur"), c("Employé", "Chômeur")))
#' chi_deux(dataE)
#'
#' # Exemple avec des vecteurs de fréquences observées et théoriques
#' obs <- c(50, 100, 120)
#' theo <- c(0.2, 0.3, 0.5)
#' chi_deux(obs, p = theo)
#' @export
chi_deux <- function (x, y = NULL, correct = TRUE, p = rep(1/length(x), length(x)),
                      rescale.p = FALSE, simulate.p.value = FALSE, B = 2000) {
  # Déclaration des noms des données
  NOM_DONNEES <- deparse(substitute(x))

  # Conversion des data frames en matrices
  if (is.data.frame(x))
    x <- as.matrix(x)

  # Conversion des matrices 1D en vecteurs
  if (is.matrix(x)) {
    if (min(dim(x)) == 1L)
      x <- as.vector(x)
  }

  # Gestion des tableaux de contingence et des vecteurs de fréquences observées et attendues
  if (!is.matrix(x) && !is.null(y)) {
    if (length(x) != length(y))
      stop("'x' et 'y' doivent avoir la même longueur")
    NOM_DONNEES2 <- deparse(substitute(y))
    nom_x <- if (length(NOM_DONNEES) > 1L || nchar(NOM_DONNEES, "w") > 30) "" else NOM_DONNEES
    nom_y <- if (length(NOM_DONNEES2) > 1L || nchar(NOM_DONNEES2, "w") > 30) "" else NOM_DONNEES2
    OK <- complete.cases(x, y)
    x <- factor(x[OK])
    y <- factor(y[OK])
    if ((nlevels(x) < 2L) || (nlevels(y) < 2L))
      stop("'x' et 'y' doivent avoir au moins 2 niveaux")
    x <- table(x, y)
    names(dimnames(x)) <- c(nom_x, nom_y)
    NOM_DONNEES <- paste(paste(NOM_DONNEES, collapse = "\n"), "et", paste(NOM_DONNEES2, collapse = "\n"))
  }

  # Vérification des valeurs négatives et des valeurs manquantes
  if (any(x < 0) || anyNA(x))
    stop("toutes les entrées de 'x' doivent être non négatives et finies")
  if ((n <- sum(x)) == 0)
    stop("au moins une entrée de 'x' doit être positive")

  # Configuration pour la simulation de la valeur p
  if (simulate.p.value) {
    setMETH <- function() METHODE <<- paste(METHODE, "avec valeur p simulée\n\t (basée sur", B, "réplicats)")
    presque.1 <- 1 - 64 * .Machine$double.eps
  }

  # Calcul des fréquences attendues et de la statistique du khi-deux pour les tableaux de contingence
  if (is.matrix(x)) {
    METHODE <- "Test du khi-deux de Pearson"
    nr <- as.integer(nrow(x))
    nc <- as.integer(ncol(x))
    if (is.na(nr) || is.na(nc) || is.na(nr * nc))
      stop("nrow(x) ou ncol(x) invalide", domain = NA)
    sr <- rowSums(x)
    sc <- colSums(x)
    E <- outer(sr, sc) / n
    v <- function(r, c, n) c * r * (n - r) * (n - c) / n^3
    V <- outer(sr, sc, v, n)
    dimnames(E) <- dimnames(x)

    # Simulation de la valeur p si demandé
    if (simulate.p.value && all(sr > 0) && all(sc > 0)) {
      setMETH()
      tmp <- .Call(C_chisq_sim, sr, sc, B, E)
      STATISTIQUE <- sum(sort((x - E)^2 / E, decreasing = TRUE))
      PARAMETRE <- NA
      PVAL <- (1 + sum(tmp >= presque.1 * STATISTIQUE)) / (B + 1)
    } else {
      if (simulate.p.value)
        warning("impossible de calculer la valeur p simulée avec des marges nulles")
      if (correct && nrow(x) == 2L && ncol(x) == 2L) {
        YATES <- min(0.5, abs(x - E))
        if (YATES > 0)
          METHODE <- paste(METHODE, "avec correction de continuité de Yates")
      } else YATES <- 0
      STATISTIQUE <- sum((abs(x - E) - YATES)^2 / E)
      PARAMETRE <- (nr - 1L) * (nc - 1L)
      PVAL <- pchisq(STATISTIQUE, PARAMETRE, lower.tail = FALSE)
    }
  } else {
    # Calcul des fréquences attendues et de la statistique du khi-deux pour les vecteurs de fréquences
    if (length(dim(x)) > 2L)
      stop("'x' invalide")
    if (length(x) == 1L)
      stop("'x' doit avoir au moins 2 éléments")
    if (length(x) != length(p))
      stop("'x' et 'p' doivent avoir le même nombre d'éléments")
    if (any(p < 0))
      stop("les probabilités doivent être non négatives.")
    if (abs(sum(p) - 1) > sqrt(.Machine$double.eps)) {
      if (rescale.p)
        p <- p / sum(p)
      else stop("les probabilités doivent totaliser 1.")
    }
    METHODE <- "Test du khi-deux pour probabilités données"
    E <- n * p
    V <- n * p * (1 - p)
    STATISTIQUE <- sum((x - E)^2 / E)
    names(E) <- names(x)

    # Simulation de la valeur p si demandé
    if (simulate.p.value) {
      setMETH()
      nx <- length(x)
      sm <- matrix(sample.int(nx, B * n, TRUE, prob = p), nrow = n)
      ss <- apply(sm, 2L, function(x, E, k) {
        sum((table(factor(x, levels = 1L:k)) - E)^2 / E)
      }, E = E, k = nx)
      PARAMETRE <- NA
      PVAL <- (1 + sum(ss >= presque.1 * STATISTIQUE)) / (B + 1)
    } else {
      PARAMETRE <- length(x) - 1
      PVAL <- pchisq(STATISTIQUE, PARAMETRE, lower.tail = FALSE)
    }
  }

  # Attribution des noms et retour des résultats
  names(STATISTIQUE) <- "Chi_deux"
  names(PARAMETRE) <- "df"
  if (any(E < 5) && is.finite(PARAMETRE))
    warning("L'approximation du khi-deux peut être incorrecte")
  structure(list(statistic = STATISTIQUE, parameter = PARAMETRE,
                 p.value = PVAL, method = METHODE, data.name = NOM_DONNEES, observed = x,
                 expected = E, residuals = (x - E) / sqrt(E), stdres = (x - E) / sqrt(V)), class = "htest")
}
