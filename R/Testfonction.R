# Création du tableau de contingence
education <- c("Primaire", "Secondaire", "Supérieur")
emploi <- c("Employé", "Chômeur")
dataE <- matrix(c(50, 100, 120, 80, 130, 20), nrow = 3, byrow = TRUE,
               dimnames = list(education, emploi))

# Conversion du tableau en data frame
df <- as.data.frame(as.table(dataE))
colnames(df) <- c("Education", "Emploi", "Nombre")

# Création du graphique avec ggplot2 et des couleurs plus douces
ggplot(df, aes(x = Education, y = Nombre, fill = Emploi)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Répartition de l'emploi en fonction du niveau d'éducation",
       x = "Niveau d'éducation", y = "Nombre de personnes") +
  theme_minimal() +
  scale_fill_manual(values = c("Employé" = "#32CD32", "Chômeur" = "#FF4500"))


# Application de la fonction au tableau de contingence
resultat <- chi_deux(dataE)
print(resultat)
View(resultat)
