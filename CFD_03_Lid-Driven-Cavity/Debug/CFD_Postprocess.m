%%
n=81;
m_grid=xlsread('100_CDS_81_velocity.csv', 'A1:CC81');
u_grid=xlsread('100_CDS_81_u.csv','A1:CC81');
v_grid=xlsread('100_CDS_81_v.csv','A1:CC81');
m1_grid=xlsread('1000_CDS_81_velocity.csv', 'A1:CC81');
u1_grid=xlsread('1000_CDS_81_u.csv','A1:CC81');
v1_grid=xlsread('1000_CDS_81_v.csv','A1:CC81');
m2_grid=xlsread('5000_CDS_81_velocity.csv', 'A1:CC81');
u2_grid=xlsread('5000_CDS_81_u.csv','A1:CC81');
v2_grid=xlsread('5000_CDS_81_v.csv','A1:CC81');
u3_grid=xlsread('5000_QUICK_81_u.csv','A1:CC81');
v3_grid=xlsread('5000_QUICK_81_v.csv','A1:CC81');
u4_grid=xlsread('5000_MUSCL_81_u.csv','A1:CC81');
v4_grid=xlsread('5000_MUSCL_81_v.csv','A1:CC81');
sp=52;
size=1.5;
%%
n=161;
m4_grid=xlsread('100_CDS_161_velocity.csv', 'A1:FE161');
u4_grid=xlsread('100_CDS_161_u.csv','A1:FE161');
v4_grid=xlsread('100_CDS_161_v.csv','A1:FE161');
m5_grid=xlsread('1000_CDS_161_velocity.csv', 'A1:FE161');
u5_grid=xlsread('1000_CDS_161_u.csv','A1:FE161');
v5_grid=xlsread('1000_CDS_161_v.csv','A1:FE161');
m6_grid=xlsread('5000_CDS_161_velocity.csv', 'A1:FE161');
u6_grid=xlsread('5000_CDS_161_u.csv','A1:FE161');
v6_grid=xlsread('5000_CDS_161_v.csv','A1:FE161');
u7_grid=xlsread('5000_QUICK_161_u.csv','A1:FE161');
v7_grid=xlsread('5000_QUICK_161_v.csv','A1:FE161');
u8_grid=xlsread('5000_MUSCL_161_u.csv','A1:FE161');
v8_grid=xlsread('5000_MUSCL_161_v.csv','A1:FE161');
% u5_grid=xlsread('100_QUICK_161_u.csv','A1:FE161');
% v5_grid=xlsread('100_QUICK_161_v.csv','A1:FE161');
% u6_grid=xlsread('100_MUSCL_161_u.csv','A1:FE161');
% v6_grid=xlsread('100_MUSCL_161_v.csv','A1:FE161');
sp=208;
size=1.5;
%%
figure(8)
hold on 
[X,Y] = meshgrid(1:n,1:n);
% image(m_grid);
contourf(m4_grid,'LineStyle', 'none');
colorbar
colormap('jet')
quiver(X(1:sp:end),Y(1:sp:end), u4_grid(1:sp:end), -v4_grid(1:sp:end), size, 'k','color','w');
xlim([1 n]) 
ylim([1 n])
set(gca, 'YDir','reverse')
axis square
%%
Ben_y_u_100=xlsread('0_Benchmark.csv', 'F1:G18');
Ben_x_v_100=xlsread('0_Benchmark.csv', 'H1:I18');
Ben_y_u_1000=xlsread('0_Benchmark.csv', 'P1:Q18');
Ben_x_v_1000=xlsread('0_Benchmark.csv', 'R1:S18');
Ben_y_u_5000=xlsread('0_Benchmark.csv', 'U1:V18');
Ben_x_v_5000=xlsread('0_Benchmark.csv', 'W1:X18');

% figure(25)
% plot(Ben_y_u_100(:,2),Ben_y_u_100(:,1),'o',u_grid(:,41),1:-1/80:0,'--r',u4_grid(:,81),1:-1/160:0,'--b',Ben_y_u_1000(:,2),Ben_y_u_1000(:,1),'og',u1_grid(:,41),1:-1/80:0,':r',u5_grid(:,81),1:-1/160:0,':b',Ben_y_u_5000(:,2),Ben_y_u_5000(:,1),'om',u2_grid(:,41),1:-1/80:0,'r',u6_grid(:,81),1:-1/160:0,'b');
% legend("Re=100 Benchmark","Re=100 81x81","Re=100 161x161","Re=1000 Benchmark","Re=1000 81x81","Re=1000 161x161","Re=5000 Benchmark","Re=5000 81x81","Re=5000 161x161",'Location','southeast')
% axis square
% 
% figure(26)
% plot(Ben_x_v_100(:,1),Ben_x_v_100(:,2),'o',0:1/80:1, v_grid(41,:),'--r',0:1/160:1, v4_grid(81,:),'--b',Ben_x_v_1000(:,1),Ben_x_v_1000(:,2),'og',0:1/80:1, v1_grid(41,:),':r',0:1/160:1, v5_grid(81,:),':b',Ben_x_v_5000(:,1),Ben_x_v_5000(:,2),'om',0:1/80:1, v2_grid(41,:),'r',0:1/160:1, v6_grid(81,:),'b');
% legend("Re=100 Benchmark","Re=100 81x81","Re=100 161x161","Re=1000 Benchmark","Re=1000 81x81","Re=1000 161x161","Re=5000 Benchmark","Re=5000 81x81","Re=5000 161x161");
% axis square


figure(25)
plot(Ben_y_u_5000(:,2),Ben_y_u_5000(:,1),'o',u2_grid(:,41),1:-1/80:0,u3_grid(:,41),1:-1/80:0,u4_grid(:,41),1:-1/80:0);
legend("Benchmark","CDS","QUICK","MUSCL",'Location','southeast')
axis square

figure(26)
plot(Ben_x_v_5000(:,1),Ben_x_v_5000(:,2),'o',0:1/80:1, v2_grid(41,:),0:1/80:1, v3_grid(41,:),0:1/80:1, v4_grid(41,:));
legend("Benchmark","CDS","QUICK","MUSCL");
axis square

