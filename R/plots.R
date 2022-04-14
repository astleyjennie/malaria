# # Plots ####

# Sort data ####
plot<-as_tibble(as.data.frame(output)) %>% 
  mutate(PThaif = Sf + Ef + Af + Cf + E2f + Rf + Tf,
         POutf = Outf + Testf,
         PThaiv = Sv + Ev + Av + Cv + E2v + Lv + Rv + Tv,
         POutv = Outv + Testv) %>% 
         #Inc1 = c(0, diff(CInc1)),
         #Inc2 = c(0, diff(CInc2))) %>%
  pivot_longer(names_to = "variable", cols = !1)%>% 
  mutate(model = "Thailand_Model")

# Population check ####
plot %>% 
  #filter(variable %in% c("PThaif", "POutf", "PThaiv", "POutv")) %>%
  filter(variable %in% c("PThaif")) %>% 
  group_by(variable, model) %>% #esquisse::esquisser()
  ggplot()+
  geom_line(aes(x = time, y=value, colour = variable))+
  theme_minimal() +
  labs(title = "Population check", y =("Populations"))
