# How to set a filter configuration?

# By indicating profile
createFilterSetting(name = "default")

# 'name' parameter has priority over the others
createFilterSetting(name = "default", type = ".definerFilter", radius = c(3, 5, 5))

# To avoid 'name' priority, it must be set as NULL
createFilterSetting(type = ".definerFilter", radius = c(3, 5, 5, 3), times = 1, tolerance = c(0.2, 0.1),
                    name = NULL)
