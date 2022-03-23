library(tidyverse)
get_dmg_decay_fit <- function(df, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  df_dx_fwd <- df %>%
    select(tax_name, label, starts_with("Dx+")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_fwd", c(-tax_name, -label)) %>%
    mutate(x = gsub("Dx\\+", "", type)) %>%
    select(-type)

  df_dx_rev <- df %>%
    select(tax_name, label, starts_with("Dx-")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_rev", c(-tax_name, -label)) %>%
    mutate(x = gsub("Dx\\-", "", type)) %>%
    select(-type)

  df_dx_std_fwd <- df %>%
    select(tax_name, label, starts_with("Dx_std+")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_std_fwd", c(-tax_name, -label)) %>%
    mutate(x = gsub("Dx_std\\+", "", type)) %>%
    select(-type)

  df_dx_std_rev <- df %>%
    select(tax_name, label, starts_with("Dx_std-")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_std_rev", c(-tax_name, -label)) %>%
    mutate(x = gsub("Dx_std\\-", "", type)) %>%
    select(-type)


  df_fit_fwd <- df %>%
    select(tax_name, label, starts_with("f+")) %>%
    pivot_longer(names_to = "type", values_to = "f_fwd", c(-tax_name, -label)) %>%
    mutate(x = gsub("f\\+", "", type)) %>%
    select(-type)

  df_fit_rev <- df %>%
    select(tax_name, label, starts_with("f-")) %>%
    pivot_longer(names_to = "type", values_to = "f_rev", c(-tax_name, -label)) %>%
    mutate(x = gsub("f\\-", "", type)) %>%
    select(-type)

  dat <- df_dx_fwd %>%
    inner_join(df_dx_rev) %>%
    inner_join(df_dx_std_fwd) %>%
    inner_join(df_dx_std_rev) %>%
    inner_join(df_fit_fwd) %>%
    inner_join(df_fit_rev) %>%
    inner_join(df %>% select(label) %>% distinct()) %>%
    mutate(x = as.numeric(x)) %>%
    filter(x <= pos) %>%
    inner_join(kapk_cdata %>% select(label, member_unit)) %>%
    rowwise() %>%
    mutate(
      Dx_fwd_min = Dx_fwd - Dx_std_fwd,
      Dx_fwd_max = Dx_fwd + Dx_std_fwd,
      Dx_rev_min = Dx_rev - Dx_std_rev,
      Dx_rev_max = Dx_rev + Dx_std_rev
    )

  fwd_max <- dat %>%
    group_by(as.character(x)) %>%
    summarise(val = mean(Dx_std_fwd) + sd(Dx_std_fwd)) %>%
    pull(val) %>%
    max()
  fwd_min <- dat %>%
    group_by(as.character(x)) %>%
    summarise(val = mean(Dx_std_fwd) - sd(Dx_std_fwd)) %>%
    pull(val) %>%
    min()
  rev_max <- dat %>%
    group_by(as.character(x)) %>%
    summarise(val = mean(Dx_std_rev) + sd(Dx_std_rev)) %>%
    pull(val) %>%
    max()
  rev_min <- dat %>%
    group_by(as.character(x)) %>%
    summarise(val = mean(Dx_std_rev) - sd(Dx_std_rev)) %>%
    pull(val) %>%
    min()

  y_max <- ifelse(fwd_max >= rev_max, fwd_max, rev_max)
  y_min <- ifelse(fwd_min <= rev_min, fwd_min, rev_min)
  y_max <- 0.5
  y_min <- -0.01
  cat(paste("fwd max:", fwd_max, "fwd min:", fwd_min, "\n"))
  cat(paste("rev max:", rev_max, "rev min:", rev_min, "\n"))
  if (orient == "fwd") {
    ggplot() +
      geom_ribbon(data = dat, aes(x, ymin = Dx_fwd_min, ymax = Dx_fwd_max, group = interaction(label, tax_name)), alpha = 0.1, fill = "black") +
      geom_line(data = dat, aes(x, Dx_fwd, group = interaction(label, tax_name)), color = "black") +
      # geom_point(data = dat, aes(x, f_fwd), alpha = .3, size = 2, fill = "black") +
      # stat_summary(data = dat, aes(x, f_fwd), fun.data = mean_sd, geom = "ribbon", alpha = .3, size = 1, fill = "black") +
      # stat_summary(data = dat, aes(x, f_fwd), fun.data = mean_sd, geom = "errorbar", alpha = .3, size = 1, fill = "black") +
      # stat_summary(data = dat, aes(x, f_fwd), fun = mean, geom = "line", size = 0.8, color = "black") +
      xlab("Position") +
      ylab("Frequency") +
      scale_y_continuous(limits = c(y_min, y_max), breaks = p_breaks) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  } else {
    ggplot() +
      geom_ribbon(data = dat, aes(x, ymin = Dx_rev_min, ymax = Dx_rev_max, group = interaction(label, tax_name)), alpha = 0.1, fill = "black") +
      geom_path(data = dat, aes(x, Dx_rev, group = interaction(label, tax_name)), color = "black") + # geom_point(data = dat, aes(x, f_rev), alpha = .3, size = 2, fill = "black") +
      # stat_summary(data = dat, aes(x, f_rev), fun.data = mean_sd, geom = "ribbon", alpha = .3, size = 1, fill = "black") +
      # stat_summary(data = dat, aes(x, f_rev), fun.data = mean_sd, geom = "errorbar", alpha = .3, size = 1, fill = "black") +
      # stat_summary(data = dat, aes(x, f_rev), fun = mean, geom = "line", size = 0.8, color = "black") +
      xlab("Position") +
      ylab("Frequency") +
      scale_x_continuous(trans = "reverse") +
      scale_y_continuous(limits = c(y_min, y_max), position = "right", breaks = p_breaks) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
}
