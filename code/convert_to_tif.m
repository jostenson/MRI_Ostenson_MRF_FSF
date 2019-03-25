%% convert png images to tiff of specified dpi

clear, close all, clc;

%%

data_dir = '../figures/';

dpi1 = 300;
dpi2 = 1200;

new_type = 'tif';

fn1{1} = 'Figure4.png';
fn1{2} = 'Figure7.png';
fn1{3} = 'Figure8.png';
fn1{4} = 'Figure9.png';
fn1{5} = 'Figure10.png';
fn1{6} = 'Figure11.png';

% fn2{1} = 'Figure1.png';
fn2{1} = 'Figure2.png';
fn2{2} = 'Figure3.png';
fn2{3} = 'Figure5.png';
fn2{4} = 'Figure6.png';


%%

n_1 = numel(fn1);
n_2 = numel(fn2);

for ii = 1:n_1

    my_img = imread( [data_dir fn1{ii}] );
    my_name = [fn1{ii}(1:end-3) new_type];
    imwrite( my_img, [data_dir my_name], 'Resolution',dpi1 );
    
end

for ii = 1:n_2

    my_img = imread( [data_dir fn2{ii}] );
    my_name = [fn2{ii}(1:end-3) new_type];
    imwrite( my_img, [data_dir my_name], 'Resolution',dpi2 );
    
end

