import json

npq_params = {}
la = [1,2,3,4]
npq_params['query_column_ranges'] = [la]
print(npq_params)
with open('a.json', 'w') as wfd:
    json.dump(npq_params, wfd)

print("done")
