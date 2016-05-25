
#loop through county names and pull data off site
fn = "http://projects-ca.statewater.org/search/county?value=sonoma"

import pandas

data = pandas.read_json(fn)

