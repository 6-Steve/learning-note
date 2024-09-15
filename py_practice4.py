#判断这一天是这一年的第几天

#获取输入
# year=int(input("请输入年份："))
# month=int(input("请输入月份："))
# day=int(input("请输入日期："))

print("请输入年份：")
year=int(input())

res=0
#创建list存储每个月份的天数（正常年的）
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

#计算天数   0~2
for i in range(month-1):
    res+=days_in_month[i]
#判断闰年
if (year % 4 == 0 and year % 100 != 0) or year % 400 == 0:
    if month > 2:
        res += 1

res=res+day

# print("%d-%d-%d是这一年的第%d天"%(year,month,day,res))
print(f"{year}-{month}-{day}是这一年的第{res}天")
